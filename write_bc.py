#rewrite of write_pc.pro in python
from read_mesa import read_mesa
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from utils import smooth, find_nearest

def find_bc(profile, mc_input):
    s    = read_mesa(profile)
    r    = np.flipud(s.radius)
    m    = np.flipud(s.mass)
    rho  = np.flipud(s.rho)
    p    = np.flipud(s.pressure)

    drhodr = np.gradient(rho,r)
    drhodr_smooth = smooth(drhodr)

    dpdr = np.gradient(p,r)
    dpdr_smooth = smooth(dpdr)

    ic = find_nearest(m, mc_input)
    bc = {'r'      : r[ic],
          'm'      : m[ic],
          'rho'    : rho[ic],
          'p'      : p[ic],
          'drhodr' : drhodr[ic],
          'dpdr'   : dpdr[ic]}

    return bc




