#python rewrite of p_profile_ambient_psmooth_extended.pro
import numpy as np
import matplotlib.pyplot as plt
from read_mesa import read_mesa
from utils import smooth, find_nearest, phantom, ohlmann
from scipy.integrate import cumtrapz, RK45


#Calculate cutoff radius:
#****IMPORTANT NOTE: MUST SET THIS MANUALLY TO BE EQUAL TO THE SAME VALUE AS THAT USED BY THE MODIFIED LANE-EMDEN SOLVER****
cutoff_by_percentage = False	#If set to True then calculate cutoff radius as a percentage of total radius using percent_try
percent_try = 5.		#approx cutoff radius in percent (used if Percentage_cutoff=1)
rc_try = 6.			#approx cutoff radius in units Rsun (used if Percentage_cutoff=0)

kernel = phantom                #Set kernel from list of kernels in utils.py

f_nt    = './out/nt.out'
f_mle   = './out/mle_profile.out'
f_param = './out/param.out'

profile = 'P12_MESA_Profile.data'

#Physical constants:
cm_Rsun = 6.957e10
g_Msun  = 1.99e33
cgs_G   = 6.6743e-8

#Read data from MLE
data_mle = np.loadtxt(f_mle).T

xi        = data_mle[0]
theta     = data_mle[1]
dthetadxi = data_mle[2]
chi       = data_mle[3]
dchidu    = data_mle[4]

#Read parameters for MLE
data_param = np.loadtxt(f_param)

dt      = data_param[0]
nnn     = data_param[1]
xi_h    = data_param[2]
xi_last = data_param[3]
rhobar  = data_param[4]
rho0    = data_param[5]
alpha   = data_param[6]
mc      = data_param[7]


s = read_mesa(profile)	#read in MESA data
r    = np.flip(s.radius,0)
m    = np.flip(s.mass,0)
rho  = 10.**np.flip(s.logRho,0)
p    = np.flip(s.pressure,0)
eint = np.flip(s.energy,0)

alpha = alpha*cm_Rsun
r = r*cm_Rsun
m = m*g_Msun
mc = mc*g_Msun

if cutoff_by_percentage:
    rc_try = 0.01 * percent_try * r[-1]    #determine cutoff radius using cutoff percentage

rc_try = rc_try * cm_Rsun

h = rc_try / kernel.radius
kernel.set_hsoft(h)

ic = find_nearest(r, rc_try)

rc   = r[ic]
mh   = m[ic] #mass m(r) at transition (cutoff) radius, different from mc which is mass of central point particle
rhoc = rho[ic]
pc   = p[ic]

r_MLE = alpha * xi



### Density profile ###

rho_MLE   = rho0 * np.power(theta,nnn)
rho_ratio = rho0 * theta[-1]**nnn / rhoc	#ratio of rho_MLE to rho_MESA at r=h

fig, axes = plt.subplots(2,3,figsize=[14,9])
axes[0,0].plot(r, rho)
axes[0,0].plot(r_MLE, rho_MLE)
axes[0,0].axvline(rc, color = 'k', dashes = [3,5])
axes[0,0].axhline(rhoc, color = 'k', dashes = [3,5])
axes[0,0].set_xscale('log')
axes[0,0].set_yscale('log')
axes[0,0].set_xlabel('log$_{10}(r)$')
axes[0,0].set_ylabel('log$_{10}(\\rho)$')

rprofile    = np.concatenate((r_MLE, r[ic:]))
rhoprofile  = np.concatenate((rho_MLE, rho[ic:]))


### Mass profile ###

# Integrate dm/dr equation from r=0 to get m(r) (serves as a check of the direct MESA output)
mmm  = [0.]
mmm2 = [0.]

for i, ri in enumerate(r[1:],1):
    dr = ri - r[i-1]
    mmm.append(mmm[i-1] + 4. * np.pi * dr * (ri**2 * rho[i] + r[i-1]**2 * rho[i-1]) / 2.)
    mmm2.append(mmm2[i-1] + 4. * np.pi * dr * r[i-1]**2 * rho[i-1])

mmmh=mmm[ic]    # mass at transition radius


# Integrate dm/dr equation for polytrope from r=0
mMLE = [0.]
for i, xii in enumerate(xi[1:], 1):
    dxi= xii - xi[i-1]
    mMLE.append(mMLE[i-1] + 4. * np.pi * dxi * alpha**3 * rho0 * (theta[i]**nnn * xii**2 + theta[i-1]**nnn * xi[i-1]**2) / 2.)

#mprofile = [0.]
#for i, ri in enumerate(rprofile[1:], 1):
#    dr = ri - rprofile[i-1]
#    mprofile.append(mprofile[i-1] + 4. * np.pi * dr * (rhoprofile[i] * ri**2 + rhoprofile[i-1] * rprofile[i-1]**2) / 2.)

#mprofile = cumtrapz(4. * np.pi * rhoprofile * rprofile**2, rprofile, initial=0.)

# Now integrate dm/dr equation from MLE solution from r=h with initial condition m(h) = m_c = m_MESA(h) to get m(r) profile
mnew = [mh]
for i, xii in reversed(list(enumerate(xi))[:-1]):
    dxi = xii - xi[i+1]
    mnew.insert(0, mnew[0] + dxi * 4. * np.pi * alpha**3 * rho0 * (theta[i]**nnn * xii**2 + theta[i+1]**nnn * xi[i+1]**2) / 2.)


mMLE = np.array(mMLE)
mnew = np.array(mnew)

print mc/g_Msun, (mc + mMLE[-1])/g_Msun, mh/g_Msun, (mc + mMLE[-1]) / mh 

axes[0,2].plot(r, m)
axes[0,2].plot(r_MLE, mc + mMLE)
#axes[0,2].plot(r_MLE, mnew)
axes[0,2].axvline(rc, color = 'k', dashes = [3,5])
axes[0,2].axhline(mh, color = 'k', dashes = [3,5])
axes[0,2].set_xscale('log')
axes[0,2].set_xlabel('log$_{10}(r)$')
axes[0,2].set_ylabel('M$_\\mathrm{internal}$')

### drho/dr profile ###

drhodr        = np.gradient(rho,r)
drhodr_smooth = smooth(drhodr)

drhodrc        = drhodr[ic]
drhodr_smoothc = drhodr_smooth[ic]

drhodr_MLE          = nnn / alpha * rho0 * theta**(nnn-1.) * dthetadxi
drhodr_ratio        = nnn / alpha * rho0 * theta[-1]**(nnn-1.) * dthetadxi[-1] / drhodrc        # ratio of drhodr_MLE to drhodr at cut radius
drhodr_ratio_smooth = nnn / alpha * rho0 * theta[-1]**(nnn-1.) * dthetadxi[-1] / drhodr_smoothc     # ratio of drhodr_MLE to drhodr_smooth at cut radius

drhodr_MLE_renorm   = drhodr_MLE / drhodr_ratio_smooth #renormalize p_MLE to equal p at r=h

axes[0,1].plot(r, -drhodr)
axes[0,1].plot(r, -drhodr_smooth[:-1])
axes[0,1].plot(r_MLE, -drhodr_MLE_renorm)
axes[0,1].axvline(rc, color = 'k', dashes = [3,5])
axes[0,1].axhline(-drhodr_smoothc, color = 'k', dashes = [3,5])
axes[0,1].set_xscale('log')
axes[0,1].set_yscale('log')
axes[0,1].set_xlabel('log$_{10}(r)$')
axes[0,1].set_ylabel('log$_{10}(d\\rho/dr)$')

### pressure profile ###

# Attempt to integrate dp/dr equation to obtain p
ppp = [pc]
for i, ri in reversed(list(enumerate(r))[1:ic+1]):# integrate backwards in radius from r=h
    dr = ri - r[i-1]
    ppp.insert(0, ppp[0] + cgs_G * m[i] * rho[i] / ri**2 * dr)


# pressure from MLE solution
p_MLE = alpha**2 * 4. * np.pi * cgs_G / (nnn + 1.) * rho0**2 * theta**(nnn + 1.)

# Now integrate dp/dr equation using MLE solution from r=h with initial condition p(h) = p_c = p_MESA(h) to get p(r) profile
pnew = [pc]
mc_mod = mc - mMLE[-1]#/g_Msun # subtract off extra mass from MLE profile for self-consistency

for i, xii in reversed(list(enumerate(xi))[:-1]):
    uuu, gc, gcdu = kernel.kernel(alpha * xii)
    gc = gc * cgs_G * mc_mod * alpha * xii / kernel.h**3
    dxi = xii - xi[i+1]
    pnew.insert(0, pnew[0] - dxi * alpha * cgs_G * mMLE[i] * rho0 * theta[i]**nnn / (xii**2 * alpha**2) - dxi * alpha * gc * rho0 * theta[i]**nnn)

p_MLE_renorm = p_MLE / (p_MLE[-1] / pc)

axes[1,0].plot(r, p)
#axes[1,1].plot(r[:ic+1], ppp)
axes[1,0].plot(r_MLE, p_MLE_renorm)
axes[1,0].plot(r_MLE, pnew)
axes[1,0].axvline(rc, color = 'k', dashes = [3,5])
axes[1,0].axhline(pc, color = 'k', dashes = [3,5])
axes[1,0].set_xscale('log')
axes[1,0].set_yscale('log')
axes[1,0].set_xlabel('log$_{10}(r)$')
axes[1,0].set_ylabel('log$_{10}(P)$')

### dpdr profile ###

dpdr = np.gradient(p,r)
dpdr_smooth = smooth(dpdr)
dpdrc= dpdr[ic]
dpdr_smoothc= dpdr_smooth[ic]

#obtain dpdr directly from modified Lane-Emden (MLE) solution by analytical differentiation of pressure
dpdr_MLE = alpha * 4. * np.pi * cgs_G * rho0**2 * theta**nnn *dthetadxi		#using modified Lane-Emden solution

#obtain dpdr from integrating dm/dr equation to get mMLE and then using hydrostic eqbm equation
dpdr_MLE_integrate = -cgs_G * (mMLE + mc) * rho0 * theta**nnn / (alpha**2 * xi**2)

#obtain dpdr by numerical differentiation of the MLE pressure profile
dpdr_MLE_deriv = np.gradient(pnew, r_MLE)

#Now obtain dp/dr using MLE soln shifted so that m(h) = m_c = m_MESA(h)
dpdrnew = []
for i, xii in reversed(list(enumerate(xi))):
    uuu, gc, gcdu = kernel.kernel(alpha * xii)
    gc = gc * cgs_G * mc_mod * alpha * xii / kernel.h**3
    dpdrnew.insert(0, -(cgs_G*mMLE[i] / (xii**2 * alpha**2) + gc) * rho0 * theta[i]**nnn)

dpdrnew = np.array(dpdrnew)

axes[1,1].plot(r, -dpdr)
axes[1,1].plot(r, -dpdr_smooth[:-1])
axes[1,1].plot(r_MLE, -dpdrnew)
axes[1,1].axvline(rc, color = 'k', dashes = [3,5])
axes[1,1].axhline(dpdr_smoothc, color = 'k', dashes = [3,5])
axes[1,1].set_xscale('log')
axes[1,1].set_yscale('log')
axes[1,1].set_xlabel('log$_{10}(r)$')
axes[1,1].set_ylabel('log$_{10}(dP/dr)$')

rprofile    = np.flip(np.concatenate((r_MLE, r[ic:])),0)
rhoprofile  = np.flip(np.concatenate((rho_MLE, rho[ic:])),0)
mprofile    = np.flip(np.concatenate((mnew, m[ic:])) - mc,0)
presprofile = np.flip(np.concatenate((pnew, p[ic:])),0)
tempprofile = np.flip(np.zeros(len(rprofile)),0)
eniprofile  = presprofile / (0.67 * rhoprofile)

### E_int profile ###
#Ideal equation of state, will not match input stellar profile

axes[1,2].plot(r, eint)
axes[1,2].plot(rprofile, eniprofile)
axes[1,2].axvline(rc, color = 'k', dashes = [3,5])
axes[1,2].set_xscale('log')
axes[1,2].set_yscale('log')
axes[1,2].set_xlabel('log$_{10}(r)$')
axes[1,2].set_ylabel('log$_{10}(E_{int})$')

plt.show()



with open('P12_Phantom_Profile_Cored.data', 'w+') as open_file:
    open_file.write('# This profile has been constructed to have a point mass core of m_c = {:5f} Msun, and a softening length of h_soft = {} Rsun.\n'.format(mc/g_Msun, kernel.h/cm_Rsun))
    open_file.write('# [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]  [{:^12}]\n'.format('Mass', 'Pressure', 'Temperature', 'Radius', 'Density', 'E_int'))
    for m, p, t, r, d,e  in zip(mprofile, presprofile, tempprofile, rprofile, rhoprofile, eniprofile):
        open_file.write('  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}\n'.format(m, p, t, r, d, e))


