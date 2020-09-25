import numpy as np
import matplotlib.pyplot as plt
import sys
from types import SimpleNamespace
from scipy.integrate import cumtrapz, simps
from scipy import optimize

from mle.read_mesa import read_mesa, read_mesa_detail
from mle.input_file_utils import read_input_file
from mle.utils import smooth, find_nearest, cm_Rsun, g_Msun, cgs_G, cgs_kb, cgs_amu
from mle.boundary_conditions import find_bc
from mle.query_mesa_eos import query_mesa_eos_solveT

# Read in the profile filenames from the command line
profile_file = sys.argv[1]
output_profile_file = sys.argv[2]

# Determine the input and output files from the MLE run
mle_output_file = profile_file.split('.')[0] + '.out'
mle_input_file = profile_file.split('.')[0] + '.in'

# Read inputs and data from MLE run
xi, theta, dthetadxi, chi, dchidu = np.loadtxt(mle_output_file).T
mle_params, _, cut_radius, mesa_eos_params = read_input_file(mle_input_file)
BCs = find_bc(profile_file, cut_radius)

# Set kernel up correctly
mle_params.kernel.set_hsoft(cut_radius / mle_params.kernel.radius)

# Convert solar units to cgs
mle_params.alpha *= cm_Rsun
mle_params.mc    *= g_Msun
mle_params.kernel.h *= cm_Rsun
BCs.r            *= cm_Rsun
BCs.m            *= g_Msun
BCs.drhodr       /= cm_Rsun
BCs.dpdr         /= cm_Rsun

# Read in MESA profile data
s = read_mesa_detail(profile_file)

# Convert solar units to cgs
s.radius *= cm_Rsun
s.mass   *= g_Msun

# Get gradients of rho and P
s.drhodr        = np.gradient(s.rho,s.radius)
drhodr_smooth   = smooth(s.drhodr)

s.dpdr         = np.gradient(s.pressure,s.radius)
dpdr_smooth    = smooth(s.dpdr)

# Reformat MESA profile quantities to be fed to the MESA EoS function
mass_fractions = np.empty((0,len(s.radius)))
for value in [e[1:] for e in mesa_eos_params.chem_id]:
	if value in dir(s):
		mass_fractions = np.vstack((mass_fractions, eval('s.'+value)))
mass_fractions = np.transpose(mass_fractions)

# Initialise cut profile - Identical to MESA profile, but with with core cut out at cut_radius
cut_profile = SimpleNamespace()
cut_profile.r = s.radius[s.radius >= BCs.r]
cut_profile.m = s.mass[s.radius >= BCs.r]
cut_profile.p = s.pressure[s.radius >= BCs.r]
cut_profile.rho = s.rho[s.radius >= BCs.r]
if 'energy' in dir(s):
	cut_profile.eint = s.energy[s.radius >= BCs.r]
	cut_profile.temp = s.temperature[s.radius >= BCs.r]
cut_profile.drhodr = s.drhodr[s.radius >= BCs.r]
cut_profile.dpdr = s.dpdr[s.radius >= BCs.r]
cut_profile.metalfrac = s.z_mass_fraction_metals[s.radius >= BCs.r]
cut_profile.hydrogenfrac = s.x_mass_fraction_H[s.radius >= BCs.r]
cut_profile.metalweig = s.abar[s.radius >= BCs.r]
cut_profile.metalchar = s.zbar[s.radius >= BCs.r]
cut_profile.massfrac = mass_fractions[s.radius >= BCs.r]

# Initialise the modified Lane-Emden profile that fills in the core
mle_profile = SimpleNamespace()
mle_profile.r = mle_params.alpha * xi
mle_profile.metalfrac = np.full((len(mle_profile.r)), cut_profile.metalfrac[0])
mle_profile.hydrogenfrac = np.full((len(mle_profile.r)), cut_profile.hydrogenfrac[0])
mle_profile.metalweig = np.full((len(mle_profile.r)), cut_profile.metalweig[0])
mle_profile.metalchar = np.full((len(mle_profile.r)), cut_profile.metalchar[0])
mle_profile.massfrac = np.full((len(mle_profile.r),mesa_eos_params.species), cut_profile.massfrac[0])

# Initialise the reconstructed profile, concatenating the MLE profile with the cut MESA profile
full_profile = SimpleNamespace()
full_profile.r = np.concatenate((mle_profile.r, cut_profile.r))

### Get modified Lane-Emden density profile ###
def MLE_density(rho0, theta, n):
    mle_rho = rho0 * np.power(theta, n)
    return mle_rho

### Get modified Lane-Emden mass profile ###
def MLE_mass(mle_r, mle_rho):
    mle_m = cumtrapz(4. * np.pi * mle_rho * mle_r**2, mle_r, initial=0.)
    return mle_m

### Get gradient of modified Lane-Emden density profile ###
def MLE_drhodr(mle_rho, alpha, n, theta, dthetadxi):
    mle_drhodr = n * mle_rho * dthetadxi / (alpha * theta)
    return mle_drhodr

### Integrate the pressure inward from the cut radius using the equations of hydrostatic equilibrium ###
def MLE_pressure(mle_r, mle_rho, mle_mass, mc, cut_radius_pres):
    # pressure from MLE solution
    chi = np.array([mle_params.kernel.kernel(r)[1] for r in mle_r])
    gc = chi * cgs_G * mc * mle_r / mle_params.kernel.h**3
    pressure = -cgs_G * mle_mass * mle_rho / mle_r**2 - mle_rho * gc
    mle_p = np.flip(cumtrapz(np.flip(pressure), np.flip(mle_r), initial=0.))
    return mle_p
		
### Get gradient of modified Lane-Emden pressure profile ###
def MLE_dpdr(mle_r, mle_p):
    # This is the analytical form, but it doesn't match up for some reason
    # dpdr = 4. * np.pi * cgs_G * alpha * rho0 * dthetadxi * mle_rho
    mle_dpdr = np.gradient(mle_p, mle_r)
    return mle_dpdr

### Calculate specific internal energy and temperature from the MESA EoS ###
def MESA_eos(metalfrac, hydrogenfrac, metalweig, metalchar, species, chem_id, net_iso, massfrac, density, pressure, temperature_guess, msg):
	#arrays for results
	eint = np.array([])
	temp = np.array([])

	#flips all the necessary arrays to loop from outside to inside
	metalfrac_flip = np.flip(metalfrac)
	hydrogenfrac_flip = np.flip(hydrogenfrac)
	metalweig_flip = np.flip(metalweig)
	metalchar_flip = np.flip(metalchar)
	massfrac_flip = np.flipud(massfrac)
	density_flip = np.flip(density)
	pressure_flip = np.flip(pressure)

	#cycle on stellar profile 
	for i, value in enumerate(density_flip):
		#specific internal energy and temperature
		eint_sol, pres_sol, temp_sol = query_mesa_eos_solveT(metalfrac_flip[i], hydrogenfrac_flip[i], metalweig_flip[i], metalchar_flip[i], species, chem_id, net_iso, massfrac_flip[i], np.log10(density_flip[i]), np.log10(pressure_flip[i]), np.log10(temperature_guess))

		#print to screen
		message = [[msg+' profile shell num', 'Specific internal energy', 'Pressure', 'Temperature'],\
                   [str(len(pressure) - i)+'/'+str(len(pressure)), str(eint_sol), str(pres_sol), str(temp_sol)]]
		col_width = max(len(word) for row in message for word in row) + 2
		for row in message:
			print("".join(word.ljust(col_width) for word in row))
		
		#updates temperature initial values
		temperature_guess = temp_sol

		#save to array
		eint = np.append(eint, eint_sol)
		temp = np.append(temp, temp_sol)

	#flips the results so that they are from the inside to the outside
	eint = np.flip(eint)
	temp = np.flip(temp)

	return eint, temp

### Calculate specific internal energy and temperature from the ideal EoS ###
def ideal_eos(rho, pres, X, Z):
    mu = (2. * X + 0.75 * (1-X-Z) + 0.56*Z)**(-1)
    eint = pres / (rho * 0.67)
    temp = pres * mu * cgs_amu / (rho * cgs_kb)

    return eint, temp

### Create profiles ###
#general quantities
mle_profile.rho = MLE_density(mle_params.rho0, theta, mle_params.n)
full_profile.rho = np.concatenate((mle_profile.rho, cut_profile.rho))

mle_profile.m = MLE_mass(mle_profile.r, mle_profile.rho)
full_profile.m = np.concatenate((mle_profile.m + mle_params.mc, cut_profile.m))

mle_profile.drhodr = MLE_drhodr(mle_profile.rho, mle_params.alpha, mle_params.n, theta, dthetadxi)
full_profile.drhodr = np.concatenate((mle_profile.drhodr, cut_profile.drhodr))

mle_profile.p = MLE_pressure(mle_profile.r, mle_profile.rho, mle_profile.m, mle_params.mc, cut_profile.p[0])
full_profile.p = np.concatenate((mle_profile.p + BCs.p, cut_profile.p))

mle_profile.dpdr = MLE_dpdr(mle_profile.r, mle_profile.p)
full_profile.dpdr = np.concatenate((mle_profile.dpdr, cut_profile.dpdr))

#mesa from interpolated T
mle_profile_temp_init = s.temperature[0]
mle_profile.eint, mle_profile.temp = MESA_eos(mle_profile.metalfrac, mle_profile.hydrogenfrac, mle_profile.metalweig, mle_profile.metalchar, mesa_eos_params.species, mesa_eos_params.chem_id, mesa_eos_params.net_iso, mle_profile.massfrac, mle_profile.rho, mle_profile.p + BCs.p, mle_profile_temp_init, 'MLE')
 #no need to recalculate the external profile, but keep this just for testing
#cut_profile_temp_init = s.temperature[-1]
#cut_profile.eint, cut_profile.temp = MESA_eos(cut_profile.metalfrac, cut_profile.hydrogenfrac, cut_profile.metalweig, cut_profile.metalchar, mesa_eos_params.species, mesa_eos_params.chem_id, mesa_eos_params.net_iso, cut_profile.massfrac, cut_profile.rho, cut_profile.p, cut_profile_temp_init, 'Outer')
full_profile.eint = np.concatenate((mle_profile.eint, cut_profile.eint))
full_profile.temp = np.concatenate((mle_profile.temp, cut_profile.temp))

#ideal
 #no need to calculate the ideal profile, but keep this just for testing
mle_profile.eint_ideal, mle_profile.temp_ideal = ideal_eos(mle_profile.rho, mle_profile.p + BCs.p, mle_profile.hydrogenfrac, mle_profile.metalfrac)
cut_profile.eint_ideal, cut_profile.temp_ideal = ideal_eos(cut_profile.rho, cut_profile.p, cut_profile.hydrogenfrac, cut_profile.metalfrac)
full_profile.eint_ideal = np.concatenate((mle_profile.eint_ideal, cut_profile.eint_ideal))
full_profile.temp_ideal = np.concatenate((mle_profile.temp_ideal, cut_profile.temp_ideal))

### Create plot ###
fig, axes = plt.subplots(2,3,figsize=[14,9])
fig.tight_layout(pad=4.0)
axes = axes.flatten()

ylabels = ['log$_{10}(P)$', 'log$_{10}(E_{int})$', 'log$_{10}(\\rho)$', 'M$_\\mathrm{internal}$', 'log$_{10}(d\\rho/dr)$', 'log$_{10}(dP/dr)$']

for i, ax in enumerate(axes[:]):
    ax.set_xlabel('log$_{10}(r)$')
    ax.set_ylabel(ylabels[i])
    ax.set_xscale('log')
    ax.set_yscale('log') if i != 3 else ax.set_yscale('linear')
    ax.axvline(BCs.r, color = 'k', dashes = [3,5])

axes[0].plot(s.radius, s.pressure, linestyle = '-')
axes[0].plot(full_profile.r, full_profile.p, linestyle = '--')

axes[1].plot(s.radius, s.energy, linestyle = '-')
axes[1].plot(full_profile.r, full_profile.eint, linestyle = '--')

axes[2].plot(s.radius, s.rho, linestyle = '-')
axes[2].plot(full_profile.r, full_profile.rho, linestyle = '--')

axes[3].plot(s.radius, s.mass, linestyle = '-')
axes[3].plot(full_profile.r, full_profile.m, linestyle = '--')
axes[3].axhline(BCs.m, color = 'k', dashes = [3,5])

axes[4].plot(s.radius, -s.drhodr, linestyle = '-')
axes[4].plot(full_profile.r, -full_profile.drhodr, linestyle = '--')
axes[4].axhline(-BCs.drhodr, color = 'k', dashes = [3,5])

axes[5].plot(s.radius, -s.dpdr, linestyle = '-')
axes[5].plot(full_profile.r, -full_profile.dpdr, linestyle = '--')
axes[5].axhline(-BCs.dpdr, color = 'k', dashes = [3,5])

plt.show()

### Save file ###
header = "This profile has been constructed to have a point mass core of m_c = {:5.3f} Msun, and a softening length of h_soft = {} Rsun.\n".format(mle_params.mc/g_Msun, mle_params.kernel.h/cm_Rsun)

titles_fmt = "[{:^12}]  " * 6
titles = titles_fmt[:-2].format('Mass', 'Pressure', 'Temperature', 'Radius', 'Density', 'E_int')

data_fmt = "  %10.8E"
data_array = np.array([
    full_profile.m - mle_params.mc,
    full_profile.p,
    full_profile.temp,
    full_profile.r,
    full_profile.rho,
    full_profile.eint
]).T

np.savetxt(output_profile_file, np.flipud(data_array), fmt=data_fmt, header=header+titles, delimiter="")
