from read_mesa import read_mesa
import numpy as np
from types import SimpleNamespace
from collections import namedtuple
import matplotlib.pyplot as plt
from utils import smooth, find_nearest

def find_bc(profile, rc_input):
    s    = read_mesa(profile)
    r    = np.flipud(s.radius)
    m    = np.flipud(s.mass)
    rho  = np.flipud(s.rho)
    p    = np.flipud(s.pressure)

    drhodr = np.gradient(rho,r)
    drhodr_smooth = smooth(drhodr)

    dpdr = np.gradient(p,r)
    dpdr_smooth = smooth(dpdr)

    ic = find_nearest(r, rc_input)
    bc = SimpleNamespace(
        r = r[ic],
        m = m[ic],
        rho = rho[ic],
        p = p[ic],
        drhodr = drhodr[ic],
        dpdr = dpdr[ic]
    )

    return bc