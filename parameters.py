from utils import *

# Initialise physical constants
g_Msun  = 1.989e33
cm_Rsun = 6.957e10
cgs_G   = 6.6743e-8

# Initialise timestepping parameters
Nmax      = 100000  # max number of points
eps_mc    = 0.00001 # sets tolerance on mc
i_dt      = 2.e-4   # iteration dt
p_dt      = 2.e-5   # profile dt
tmin      = 1.e-10  # must start at finite value to avoid singularity

# Initialise physical parameters
nnn          = 3      # polytropic index
alpha_init   = 1.0   # radius parameter
rho0_init    = 0.5 # central density, units g/cm^3
mc_on_mh_min = 0.5    # lower guess for ratio of core mass to cut mass
mc_on_mh_max = 1.0    # upper guess for ratio of core mass to cut mass

# Set kernel type
kernel = phantom

# Cut conditions
cut_radius   = 1      # Units are Rsun
