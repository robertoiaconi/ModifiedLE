#python rewrite of p_profile_ambient_psmooth_extended.pro
import numpy as np
import matplotlib.pyplot as plt
from read_mesa import read_mesa
from utils import smooth, find_nearest
from scipy.integrate import cumtrapz, odeint
from parameters import *
from eos import *
import sys

profile = sys.argv[1]
f_mle = 'out/{}.out'.format(profile.split('.')[0])
f_param = 'out/{}.param'.format(profile.split('.')[0])
outprofile = sys.argv[2]

#Read data from MLE
data_mle = np.loadtxt(f_mle).T
xi, theta, dthetadxi, chi, dchidu = data_mle

#Read parameters for MLE
data_param = np.loadtxt(f_param)
dt, nnn, xi_h, xi_last, rhobar, rho0, alpha, m_cut = data_param

s = read_mesa(profile)	#read in MESA data
r    = np.flipud(s.radius)
m    = np.flipud(s.mass)
rho  = np.flipud(s.rho)
p    = np.flipud(s.pressure)
eint = np.flipud(s.energy)
t    = np.flipud(s.temperature)

alpha *= cm_Rsun
r     *= cm_Rsun
m     *= g_Msun
m_cut *= g_Msun
input_mc *= g_Msun

drhodr        = np.gradient(rho,r)
drhodr_smooth = smooth(drhodr)

dpdr         = np.gradient(p,r)
dpdr_smooth  = smooth(dpdr)

kernel.set_hsoft(input_h*cm_Rsun)

ic = find_nearest(m, m_cut)
rc   = r[ic]
mc   = m[ic] #mass m(r) at transition (cutoff) radius, different from mc which is mass of central point particle
rhoc = rho[ic]
pc   = p[ic]

X_mle    = np.flipud(s.x_mass_fraction_H)[ic]
Y_mle    = np.flipud(s.y_mass_fraction_He)[ic]
Z_mle    = np.flipud(s.z_mass_fraction_metals)[ic]

r_MLE = alpha * xi
len_MLE = len(r_MLE)
rprofile = np.concatenate((r_MLE[:-1], r[ic:]))
print(alpha, rho0)

### Density profile ###

def MLE_density(rho0, theta, nnn, rho):
    rho_MLE = rho0 * np.power(theta,nnn)
    return np.concatenate((rho_MLE[:-1], rho))

### Mass profile ###

def MLE_mass1(alpha, rho0, xi, theta, nnn, m):
    # Integrate dm/dr equation for polytrope from r=0
    mMLE = cumtrapz(4. * np.pi * alpha**2 * xi**2 * rho0 * theta**nnn, alpha * xi, initial=0.)

    # Now integrate dm/dr equation from MLE solution from r=h with initial condition m(h) = m_c = m_MESA(h) to get m(r) profile
    mnew = [mc]
    for i, xii in reversed(list(enumerate(xi))[:-1]):
        dxi = xii - xi[i+1]
        mnew.insert(0, mnew[0] + dxi * 4. * np.pi * alpha**3 * rho0 * (theta[i]**nnn * xii**2 + theta[i+1]**nnn * xi[i+1]**2) / 2.)

    print(' Core mass = {}\n MLE mass at cut = {}\n MESA mass at cut = {}\n Mass ratio at cut = {}'.format(input_mc/g_Msun, (input_mc + mMLE[-1])/g_Msun, mc/g_Msun, (input_mc + mMLE[-1]) / mc))
    return np.concatenate((mnew[:-1], m)) - input_mc, mMLE

def mass_func(m,r,rho):
    dmdr = 4 * np.pi * r**2 * rho
    return dmdr      


def MLE_mass(alpha, rho0, xi, theta, nnn, m):
    # Integrate dm/dr equation for polytrope from r=0
    mMLE = cumtrapz(4. * np.pi * alpha**2 * xi**2 * rho0 * theta**nnn, alpha * xi, initial=0.)

    r = alpha * xi * cm_Rsun
    rho = rho0 * np.power(theta,nnn)

    sol = odeint(mass_func, (0, -r[-1]), [mc], method='Radau', dense_output=True)
    print(sol)
    plt.loglog(sol.t + r[-1], sol.y[0] + input_mc)
    plt.show()
    exit()

    # Now integrate dm/dr equation from MLE solution from r=h with initial condition m(h) = m_c = m_MESA(h) to get m(r) profile
    mnew.insert(0, mnew[0] + dxi * 4. * np.pi * alpha**3 * rho0 * (theta[i]**nnn * xii**2 + theta[i+1]**nnn * xi[i+1]**2) / 2.)

    print(' Core mass = {}\n MLE mass at cut = {}\n MESA mass at cut = {}\n Mass ratio at cut = {}'.format(input_mc/g_Msun, (input_mc + mMLE[-1])/g_Msun, mc/g_Msun, (input_mc + mMLE[-1]) / mc))
    return np.concatenate((mnew[:-1], m)) - input_mc, mMLE


### drho/dr profile ###

def MLE_drhodr(alpha, rho0, xi, theta, dthetadxi, nnn, drhodr):
    drhodr_MLE          = nnn / alpha * rho0 * theta**(nnn-1.) * dthetadxi
    drhodr_ratio        = nnn / alpha * rho0 * theta[-1]**(nnn-1.) * dthetadxi[-1] / drhodr[0]        # ratio of drhodr_MLE to drhodr at cut radius

    return np.concatenate((drhodr_MLE[:-1] / drhodr_ratio, drhodr)) #renormalize p_MLE to equal p at r=h

### pressure profile ###

# Attempt to integrate dp/dr equation to obtain p
def MLE_pressure(alpha, rho0, xi, theta, dthetadxi, nnn, p, m):
    # pressure from MLE solution
    p_MLE = alpha**2 * 4. * np.pi * cgs_G / (nnn + 1.) * rho0**2 * theta**(nnn + 1.)

    # Now integrate dp/dr equation using MLE solution from r=h with initial condition p(h) = p_c = p_MESA(h) to get p(r) profile
    pnew = [pc]
    mc_mod = input_mc #- m[-1]#/g_Msun # subtract off extra mass from MLE profile for self-consistency

    for i, xii in reversed(list(enumerate(xi))[:-1]):
        uuu, gc, gcdu = kernel.kernel(alpha * xii)
        gc = gc * cgs_G * mc_mod * alpha * xii / kernel.h**3
        dxi = xii - xi[i+1]
        pnew.insert(0, pnew[0] - dxi * alpha * cgs_G * m[i] * rho0 * theta[i]**nnn / (xii**2 * alpha**2) - dxi * alpha * gc * rho0 * theta[i]**nnn)
    return np.concatenate((pnew[:-1],p))

### dpdr profile ###

def MLE_dpdr(alpha, rho0, xi, theta, dthetadxi, nnn, dpdr, m):
    #obtain dpdr directly from modified Lane-Emden (MLE) solution by analytical differentiation of pressure
    dpdr_MLE = alpha * 4. * np.pi * cgs_G * rho0**2 * theta**nnn *dthetadxi		#using modified Lane-Emden solution

    #obtain dpdr from integrating dm/dr equation to get mMLE and then using hydrostic eqbm equation
    #dpdr_MLE_integrate = -cgs_G * (m+mc) * rho0 * theta**nnn / (alpha**2 * xi**2)

    #obtain dpdr by numerical differentiation of the MLE pressure profile
    #dpdr_MLE_deriv = np.gradient(p, r)

    #Now obtain dp/dr using MLE soln shifted so that m(h) = m_c = m_MESA(h)
    dpdrnew = []
    mc_mod = input_mc #- m[-1]#/g_Msun # subtract off extra mass from MLE profile for self-consistency
    for i, xii in reversed(list(enumerate(xi))):
        uuu, gc, gcdu = kernel.kernel(alpha * xii)
        gc = gc * cgs_G * mc_mod * alpha * xii / kernel.h**3
        dpdrnew.insert(0, -(cgs_G*m[i] / (xii**2 * alpha**2) + gc) * rho0 * theta[i]**nnn)

    return np.concatenate((dpdrnew[:-1], dpdr))

def MLE_eint_temp(rho, pres, X, Z):
    mu = (2. * X + 0.75 * (1-X-Z) + 0.56*Z)**(-1)
    kb = 1.38e-16
    amu = 1.66e-24 
    eint = pres / (rho * 0.67)
    temp = pres * mu * amu / (rho * kb)

    return temp, eint

### Create profiles ###

rhoprofile = MLE_density(rho0, theta, nnn, rho[ic:])
mprofile, mMLE = MLE_mass(alpha, rho0, xi, theta, nnn, m[ic:])
drhodrprofile = MLE_drhodr(alpha, rho0, xi, theta, dthetadxi, nnn, drhodr_smooth[ic:-1])
presprofile = MLE_pressure(alpha, rho0, xi, theta, dthetadxi, nnn, p[ic:], mMLE)
dpdrprofile = MLE_dpdr(alpha, rho0, xi, theta, dthetadxi, nnn, dpdr_smooth[ic:-1], mMLE) 
tempprofile, eniprofile = query_temp(rhoprofile[:len_MLE], presprofile[:len_MLE], X_mle, Z_mle)
#tempprofile, eniprofile = MLE_eint_temp(rhoprofile[:], presprofile[:], X_mle, Z_mle)
eniprofile = np.concatenate((eniprofile[:-1], eint[ic:]))
tempprofile = np.concatenate((tempprofile[:-1], t[ic:]))
### Create plot ###

fig, axes = plt.subplots(2,3,figsize=[14,9])
axes = axes.flatten()

ylabels = ['log$_{10}(\\rho)$', 'log$_{10}(d\\rho/dr)$', 'M$_\\mathrm{internal}$', 'log$_{10}(P)$', 'log$_{10}(dP/dr)$', 'log$_{10}(E_{int})$']

for i, ax in enumerate(axes[:]):
    ax.set_xlabel('log$_{10}(r)$')
    ax.set_ylabel(ylabels[i])
    ax.set_xscale('log')
    ax.set_yscale('log') if i != 2 else ax.set_yscale('linear')
    ax.axvline(rc, color = 'k', dashes = [3,5])

axes[0].plot(r, rho)
axes[0].plot(rprofile, rhoprofile)
axes[0].axhline(rhoc, color = 'k', dashes = [3,5])

axes[1].plot(r, -drhodr)
axes[1].plot(r, -drhodr_smooth[:-1], color='tab:green')
axes[1].plot(rprofile, -drhodrprofile, color='tab:orange')
axes[1].axhline(-drhodr_smooth[ic], color = 'k', dashes = [3,5])

axes[2].plot(r, m)
axes[2].plot(rprofile, mprofile + input_mc)
axes[2].axhline(mc, color = 'k', dashes = [3,5])

axes[3].plot(r, p)
axes[3].plot(rprofile, presprofile)
axes[3].axhline(pc, color = 'k', dashes = [3,5])

axes[4].plot(r, -dpdr)
axes[4].plot(r, -dpdr_smooth[:-1], color='tab:green')
axes[4].plot(rprofile, -dpdrprofile, color='tab:orange')
axes[4].axhline(-dpdr_smooth[ic], color = 'k', dashes = [3,5])

axes[5].plot(r, eint)
axes[5].plot(rprofile, eniprofile)

plt.show()

with open(outprofile, 'w+') as open_file:
    rprofile    = np.flipud(rprofile)
    rhoprofile  = np.flipud(rhoprofile)
    mprofile    = np.flipud(mprofile)
    presprofile = np.flipud(presprofile)
    eniprofile  = np.flipud(eniprofile)
    tempprofile = np.flipud(tempprofile)

    open_file.write('# This profile has been constructed to have a point mass core of m_c = {:5f} Msun, and a softening length of h_soft = {} Rsun.\n'.format(input_mc/g_Msun, kernel.h/cm_Rsun))
    open_file.write('# [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]\n'.format('Mass', 'Pressure', 'Temperature', 'Radius', 'Density', 'E_int'))
    for m, p, t, r, d,e  in zip(mprofile, presprofile, tempprofile, rprofile, rhoprofile, eniprofile):
        open_file.write('  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}\n'.format(m, p, t, r, d, e))
