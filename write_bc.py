#rewrite of write_pc.pro in python
from read_mesa import read_mesa
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from utils import smooth, find_nearest


nsmooth_drhodr= 10	# for smoothing of drho/dr profile
nsmooth_dpdr=   10	# for smoothing of dp/dr profile

datafile= 'bc.out'
profile= 'P12_MESA_Profile.data'

cutoff_by_percentage = False	#If set to 1 then calculate cutoff radius as a percentage of total radius using percent_try

percent_try = 5.	#desired cutoff radius in percent (used if cutoff_by_percentage = True)
rc_try= 6	#desired cutoff radius in units Rsun (used if cutoff_by_percentage = False)

s    = read_mesa(profile)        #use IDL function from MESA website to read in MESA log file
r    = np.flip(s.radius,0)
m    = np.flip(s.mass,0)
rho  = 10.**np.flip(s.logRho,0)
p    = np.flip(s.pressure,0)
pgas = np.flip(s.pressure*s.pgas_div_ptotal,0) # gas pressure - can get from pgas in profile
prad = p-pgas                                  # radiation pressure

drhodr = np.gradient(rho,r)
drhodr_smooth = smooth(drhodr)

dpdr = np.gradient(p,r)
dpdr_smooth = smooth(dpdr)

if cutoff_by_percentage:
    rc_try = 0.01 * percent_try * r[-1]

#determine closest values of these quantities at the cutoff radius r_c

ic = find_nearest(r, rc_try)

rc             = r[ic]
mc             = m[ic]
rhoc           = rho[ic]
pc             = p[ic]
pgasc          = pgas[ic]
pradc          = prad[ic]
drhodrc        = drhodr[ic]
drhodr_smoothc = drhodr_smooth[ic]
dpdrc          = dpdr[ic]
dpdr_smoothc   = dpdr_smooth[ic]

with open('bc.out', 'w+') as open_file:
    open_file.write(str(rc)+' ')
    open_file.write(str(mc)+' ')
    open_file.write(str(rhoc)+' ')
    open_file.write(str(drhodr_smoothc)+' ')
    open_file.write(str(pc)+' ')
    open_file.write(str(dpdr_smoothc)+' ')




