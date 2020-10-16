import numpy as np
from types import SimpleNamespace
import matplotlib.pyplot as plt
from .read_mesa import read_mesa
from .utils import smooth, find_nearest

def find_bc(profile, rc_input):
    s = read_mesa(profile)

    drhodr = np.gradient(s.rho,s.radius)
    drhodr_smooth = smooth(drhodr)

    dpdr = np.gradient(s.pressure,s.radius)
    dpdr_smooth = smooth(dpdr)

    ic = find_nearest(s.radius, rc_input)
    bc = SimpleNamespace(
        r = s.radius[ic],
        m = s.mass[ic],
        rho = s.rho[ic],
        p = s.pressure[ic],
        drhodr = drhodr[ic],
        dpdr = dpdr[ic]
    )

    return bc
