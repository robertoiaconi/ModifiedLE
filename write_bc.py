import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline

from read_mesa import read_mesa
from utils import find_nearest
from parameters import g_Msun, cm_Rsun

# Return the boundary conditions
def find_bc(profile, rc_input):
    s    = read_mesa(profile)
    r    = np.flipud(s.radius)
    m    = np.flipud(s.mass)
    rho  = np.flipud(s.rho)
    p    = np.flipud(s.pressure)

    drhodr = np.gradient(rho,r)
    dpdr = np.gradient(p,r)

    ic = find_nearest(r, rc_input)
    BCs = namedtuple('BCs', ['r', 'm', 'rho', 'p', 'drhodr', 'dpdr'])
    bc = BCs(r      = r[ic] * cm_Rsun,
             m      = m[ic] * g_Msun,
             rho    = rho[ic],
             p      = p[ic],
             drhodr = drhodr[ic],
             dpdr   = dpdr[ic])

    return bc



