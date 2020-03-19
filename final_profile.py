import numpy as np
import matplotlib.pyplot as plt
import sys
from types import SimpleNamespace
from scipy.integrate import cumtrapz, simps

from mle.read_mesa import read_mesa
from mle.input_file_utils import read_input_file
from mle.utils import smooth, find_nearest, cm_Rsun, g_Msun, cgs_G, cgs_kb, cgs_amu
from mle.boundary_conditions import find_bc

# Read in the profile filenames from the command line
profile_file = sys.argv[1]
output_profile_file = sys.argv[2]

# Determine the input and output files from the MLE run
mle_output_file = profile_file.split('.')[0] + '.out'
mle_input_file = profile_file.split('.')[0] + '.in'

# Read inputs and data from MLE run
xi, theta, dthetadxi, chi, dchidu = np.loadtxt(mle_output_file).T
mle_params, _, cut_radius = read_input_file(mle_input_file)
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
s = read_mesa(profile_file)

# Convert solar units to cgs
s.radius *= cm_Rsun
s.mass   *= g_Msun

# Get gradients of rho and P
s.drhodr        = np.gradient(s.rho,s.radius)
drhodr_smooth   = smooth(s.drhodr)

s.dpdr         = np.gradient(s.pressure,s.radius)
dpdr_smooth    = smooth(s.dpdr)

# Initialise cut profile - Identical to MESA profile, but with with core cut out at cut_radius
cut_profile = SimpleNamespace()
cut_profile.r = s.radius[s.radius >= BCs.r]
cut_profile.m = s.mass[s.radius >= BCs.r]
cut_profile.p = s.pressure[s.radius >= BCs.r]
cut_profile.rho = s.rho[s.radius >= BCs.r]
cut_profile.drhodr = s.drhodr[s.radius >= BCs.r]
cut_profile.dpdr = s.dpdr[s.radius >= BCs.r]

# Initialise the modified Lane-Emden profile that fills in the core
mle_profile = SimpleNamespace()
mle_profile.r = mle_params.alpha * xi

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
def MLE_pressure(mle_r, mle_rho, mle_mass, mc):
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

def ideal_eos(rho, pres, X, Z):
    mu = (2. * X + 0.75 * (1-X-Z) + 0.56*Z)**(-1)
    eint = pres / (rho * 0.67)
    temp = pres * mu * cgs_amu / (rho * cgs_kb)

    return eint, temp


### Create profiles ###
mle_profile.rho = MLE_density(mle_params.rho0, theta, mle_params.n)
full_profile.rho = np.concatenate((mle_profile.rho, cut_profile.rho))

mle_profile.m = MLE_mass(mle_profile.r, mle_profile.rho)
full_profile.m = np.concatenate((mle_profile.m + mle_params.mc, cut_profile.m))

mle_profile.drhodr = MLE_drhodr(mle_profile.rho, mle_params.alpha, mle_params.n, theta, dthetadxi)
full_profile.drhodr = np.concatenate((mle_profile.drhodr, cut_profile.drhodr))

mle_profile.p = MLE_pressure(mle_profile.r, mle_profile.rho, mle_profile.m, mle_params.mc)
full_profile.p = np.concatenate((mle_profile.p + BCs.p, cut_profile.p))

mle_profile.dpdr = MLE_dpdr(mle_profile.r, mle_profile.p)
full_profile.dpdr = np.concatenate((mle_profile.dpdr, cut_profile.dpdr))

full_profile.eint, full_profile.temp = ideal_eos(full_profile.rho, full_profile.p, X=0.75, Z=0.)

### Create plot ###

fig, axes = plt.subplots(2,3,figsize=[14,9])
axes = axes.flatten()

ylabels = ['log$_{10}(\\rho)$', 'log$_{10}(d\\rho/dr)$', 'M$_\\mathrm{internal}$', 'log$_{10}(P)$', 'log$_{10}(dP/dr)$', 'log$_{10}(E_{int})$']

for i, ax in enumerate(axes[:]):
    ax.set_xlabel('log$_{10}(r)$')
    ax.set_ylabel(ylabels[i])
    ax.set_xscale('log')
    ax.set_yscale('log') if i != 2 else ax.set_yscale('linear')
    ax.axvline(BCs.r, color = 'k', dashes = [3,5])

axes[0].plot(s.radius, s.rho)
axes[0].plot(full_profile.r, full_profile.rho)
axes[0].axhline(BCs.rho, color = 'k', dashes = [3,5])

axes[1].plot(s.radius, -s.drhodr)
axes[1].plot(full_profile.r, -full_profile.drhodr)
axes[1].axhline(-BCs.drhodr, color = 'k', dashes = [3,5])

axes[2].plot(s.radius, s.mass)
axes[2].plot(full_profile.r, full_profile.m)
axes[2].axhline(BCs.m, color = 'k', dashes = [3,5])

axes[3].plot(s.radius, s.pressure)
axes[3].plot(full_profile.r, full_profile.p)
axes[3].axhline(BCs.p, color = 'k', dashes = [3,5])

axes[4].plot(s.radius, -s.dpdr)
axes[4].plot(full_profile.r, -full_profile.dpdr)
axes[4].axhline(-BCs.dpdr, color = 'k', dashes = [3,5])

axes[5].plot(s.radius, s.energy)
axes[5].plot(full_profile.r, full_profile.eint)

plt.show()


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
