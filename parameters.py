from utils import *

# Initialise physical constants
g_Msun  = 1.989e33
cm_Rsun = 6.957e10
cgs_G   = 6.6743e-8

# Initialise timestepping parameters
Nmax      = 100000  # max number of points
eps_mc    = 0.000001 # sets tolerance on mc
i_dt      = 2.e-4   # iteration dt
p_dt      = 2.e-4   # profile dt
tmin      = 1.e-10  # must start at finite value to avoid singularity

# Initialise physical parameters
nnn          = 3.     # polytropic index
alpha_init   = 10.     # radius parameter
rho0_init    = 0.005  # central density, units g/cm^3
min_mh       = 0.392  # lower guess for ratio of core mass to cut mass
max_mh       = 0.41   # upper guess for ratio of core mass to cut mass

# Set kernel type
kernel = phantom

# Cut conditions
input_mc     = 0.392        # if cutoff_by_percentage is False, units are Rsun, otherwise %
input_h      = 3

# Ivanova cut_radius = 2   : alpha = 2.3    rho0 = 0.017
# Ivanova cut_radius = 0.7 : alpha = 0.8    rho0 = 0.13
# P12     cut radius = 6.  : alpha = 5.57   rho0 = 0.00056

# Chang   cut_radius = 6.  : alpha = 6.9    rho0 = 0.00055    Mcore = 0.372
# Chang   cut_radius = 5.  : alpha = 6.11   rho0 = 0.00075    Mcore = 0.369

