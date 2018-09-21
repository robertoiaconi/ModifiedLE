import os
import sys
import numpy as np
from scipy.integrate import cumtrapz
from utils import phantom, ohlmann

# Initialise physical constants
s_Myr   = 1e6 * 365.25 * 24. * 3600.
s_Gyr   = 1e9 * 365.25 * 24. * 3600.
g_Msun  = 1.989e33
cm_Rsun = 6.957e10
cgs_G   = 6.6743e-8

# Initialise timestepping parameters
Nmax      = 100000 #max number of points
ialphamax = 100000 #max iterations to find alpha to match BCs
irho0max  = 100000 #max iterations to find rho0 to match BCs
eps_alpha = 0.001  #sets tolerance of convergence condition on alpha--smaller eps is more precise
eps2      = 0.0005 #sets fractional change in alpha after each iteration
eps_rho0  = 0.001  #sets tolerance of convergence condition on rho0--smaller eps is more precise
eps_mc    = 0.001
i_dt      = 2.e-4  #iteration dt
p_dt      = 2.e-5  #profile dt
tmin      = 1.e-10 #must start at finite value to avoid singularity

# Initialise bc file
bc_file = 'bc.out'
output_file = "out/mle_profile.out"
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

# Initialise physical parameters
nnn          = 3     #polytropic index
alpha_init   = 5.87   #radius parameter
rho0_init    = 0.0005 #central density, units g/cm^3
mc_on_mh_min = 0.9
mc_on_mh_max = 1.

# Set the desired kernel
kernel = phantom

def set_bc(bc_file, alpha):
    with open(bc_file,'r') as open_file:
        bc_out = open_file.readline().split()
	
    rcut     = float(bc_out[0])
    h        = rcut / kernel.radius
    m_h      = float(bc_out[1])
    rho_h    = float(bc_out[2])
    drhodr_h = float(bc_out[3])

    return rcut, h, m_h, rho_h, drhodr_h


def pde(t, f, alpha, xi_h, nnn, rhobar, rho0):
    xi        = t
    theta     = f[0]
    eta       = f[1]
    dthetadxi = -eta/(xi**2)

    u, chi, dchidu = kernel.kernel(alpha * xi)
    
    dfdt0 = -eta / xi**2
    dfdt1 = xi**2 * (theta**nnn + rhobar / rho0 * (chi + 1./3. * u * dchidu))

    dfdt = np.array([dfdt0, dfdt1])

    return dfdt


def rk3(t, dt, f, counter, alpha, xi_h, nnn, rhobar, rho0): # Exact copy of third order RK integrator from fortran version
    #gam1, gam2 AND gam3 ARE THE COEFFICIENTS OF THE TIMESTEPS AT WHICH dfdt IS CALCULATED
    gam1 = 8./15. 
    gam2 = 5./12.
    gam3 = 3./4.
    zet1 = -17./60.
    zet2 = -5./12.

    dfdt = pde(t, f, alpha, xi_h, nnn, rhobar, rho0)
    pdef = dfdt                     #pde FUNCTION CALCULATES VALUE OF TIME DERIVATIVE ACCORDING TO P.D.E.s
    f    = f + dt * gam1 * pdef     #FIRST GET INTERMEDIATE VALUE OF f (AT t=t+dt*gam1) USING dfdt AT t_0
    t    = t + dt * gam1            #THEN GO TO THAT TIMESTEP
    ftmp = f + dt * zet1 * pdef     #NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1) USING dfdt AT t_0
    ttmp = t + dt * zet1            #DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1)

    dfdt = pde(t, f, alpha, xi_h, nnn, rhobar, rho0)
    pdef = dfdt                     #NOW CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1
    f    = ftmp + dt * gam2 * pdef  #USE THIS TO GET ANOTHER INTERMEDIATE VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2)
    t    = ttmp + dt * gam2         #THEN GO TO THAT TIMESTEP
    ftmp = f + dt * zet2 * pdef     #NOW CALCULATE A TEMPORARY f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2) USING dfdt AT t=t_0+dt*gam1
    ttmp = t + dt * zet2            #DEFINE TEMPORARY t AS THAT TIMESTEP (t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2)

    dfdt = pde(t, f, alpha, xi_h, nnn, rhobar, rho0)
    pdef = dfdt                     #CALCULATE dfdt AT THE NEW TIMESTEP t_0+dt*gam1+dt*zet1+dt*gam2
    f    = ftmp + dt * gam3 * pdef  #USE THIS TO GET THE FINAL VALUE OF f (AT t=t_0+dt*gam1+dt*zet1+dt*gam2+dt*zet2+dt*gam3)
    t    = ttmp + dt * gam3         #THEN GO TO THAT TIMESTEP

    counter += 1                    #COUNTS THE NUMBER OF TIMES RUNGA-KUTTA ROUTINE IS EXCUTED

    return t, f, counter


def init_start(rho0, alpha):
    t       = tmin   #initial value of t
    counter = 1      #counts RK calls

    f = np.array([1., 0.]) #Initial conditions for integration

    return t, f, counter

def eval_rk(alpha, rho0, xi_h, rhobar, write=False):
    t, f, counter = init_start(rho0, alpha)
    xi, theta, eta = np.array([t]), np.array([f[0]]), np.array([f[1]])
    dt = i_dt
    if write:
        open_file = open(output_file, 'w+')
        dt = p_dt

    for i in range(Nmax):
        t, f, counter = rk3(t, dt, f, counter, alpha, xi_h, nnn, rhobar, rho0)

        xi    = np.append(xi,t)
        theta = np.append(theta,f[0])
        eta   = np.append(eta,f[1])

        if write:
            dthetadxi = -eta[-1] / xi[-1]**2
            u, chi, dchidu = kernel.kernel(alpha * xi[-1])
            open_file.write('{:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}\n'.format(xi[-1], theta[-1], dthetadxi, chi, dchidu))

        if np.isnan(f[0]):
            print 'NaNs detected, change the time step'
            exit()

        if f[0] <= 0:
            print 'reached negative theta at it = ', i
            print 'xi = ', t
            exit()
        elif xi[-1] >= xi_h:
            break
    
    if write:
        open_file.close()

    return xi, eta, theta

def rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha):
    xi, eta, theta = eval_rk(alpha, rho0, xi_h, rhobar)

    r = alpha * xi
    ircut = (np.abs(r-rcut)).argmin()
    rho = rho0 * theta[ircut]**nnn
    rhop = rho0 * theta[ircut-1]**nnn
    drhodr = (rho - rhop) / (r[ircut] - r[ircut-1])  #This could probably be improved
    return rho, drhodr

# Find how rho and drhodr vary as alpha varies around present value of alpha
def drhodalpha(rcut, rhobar, nnn, rho0, alpha, xi_h):
    h_alpha = np.sqrt(sys.float_info.epsilon) * alpha
    rhon2, drhon2 = rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha-2*h_alpha)
    rhon1, drhon1 = rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha-h_alpha)
    rhop1, drhop1 = rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha+h_alpha)
    rhop2, drhop2 = rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha+2*h_alpha)
    drho = (rhon2 - 8.*rhon1 + 8.*rhop1 - rhop2) / (12.*h_alpha)
    ddrho = (drhon2 - 8.*drhon1 + 8.*drhop1 - drhop2) / (12.*h_alpha)
    return drho, ddrho

# Find how rho and drhodr vary as rho0 varies around present value of rho0
def drhodrho0(rcut, rhobar, nnn, rho0, alpha, xi_h):
    h_rho0 = np.sqrt(sys.float_info.epsilon) * rho0
    rhon2, drhon2 = rhocut(rcut, rhobar, nnn, xi_h, rho0-2*h_rho0, alpha)
    rhon1, drhon1 = rhocut(rcut, rhobar, nnn, xi_h, rho0-h_rho0, alpha)
    rhop1, drhop1 = rhocut(rcut, rhobar, nnn, xi_h, rho0+h_rho0, alpha)
    rhop2, drhop2 = rhocut(rcut, rhobar, nnn, xi_h, rho0+2*h_rho0, alpha)
    drho = (rhon2 - 8.*rhon1 + 8.*rhop1 - rhop2) / (12.*h_rho0)
    ddrho = (drhon2 - 8.*drhon1 + 8.*drhop1 - drhop2) / (12.*h_rho0)
    return drho, ddrho

# Newton-Raphson with two variables, matching rho = rho_h and drhodr = drhodr_h
def newton_raphson(alpha, rho0, rcut, rhobar, rho_h, drhodr_h):

    # Initialise Jacobian
    J = np.array([[0.,0.],[0.,0.]])

    for i in xrange(100): # 100 is arbitrary, but if it hasn't converged by then it likely won't

        xi_h = rcut / alpha
        J[0][0], J[1][0] = drhodalpha(rcut, rhobar, nnn, rho0, alpha, xi_h) # J = [drho/dalpha    d^2rho/drdalpha]
        J[0][1], J[1][1] = drhodrho0(rcut, rhobar, nnn, rho0, alpha, xi_h)  #     [drho/drho0     d^2rho/drdrho0 ]

        J = np.linalg.inv(J) # Matrix inverse of J

        rho, drhodr = rhocut(rcut, rhobar, nnn, xi_h, rho0, alpha)     # Find rho and drho/dr at the cut radius
        alpha = alpha - np.dot(J, [rho - rho_h, drhodr - drhodr_h])[0] # Increment alpha based on gradients
        rho0 = rho0 - np.dot(J, [rho - rho_h, drhodr - drhodr_h])[1]   # Increment rho0 based on gradients

        rho_check = abs((rho_h-rho) / rho)
        drhodr_check = abs((drhodr_h - drhodr) / drhodr)

        print '\n##### Iteration {} #####'.format(i)
        print 'alpha = {}'.format(alpha)
        print 'rho0  = {}'.format(rho0)
        print 'Desired cut values : rho = {}, drhodr = {}'.format(rho_h, drhodr_h)
        print 'MLE cut values     : rho = {}, drhodr = {}'.format(rho, drhodr)
        print 'Cut rho difference    = {}'.format(rho_check)
        print 'Cut drhodr difference = {}'.format(drhodr_check)

        if (rho_check < 1e-9 and drhodr_check < 1e-9):
            break
    return alpha, rho0, xi_h


#Main function
def mle_run(kernel):
    # Initialise boundary conditions and initial guesses for alpha and rho0
    alpha = alpha_init
    rho0  = rho0_init
    rcut, h, m_h, rho_h, drhodr_h = set_bc(bc_file, alpha)

    # Set softening length of the kernel
    kernel.set_hsoft(h)

    # Guts of the code - use bisection root finder to iterate m_c until the cut mass matches m_h
    alpha, rho0, xi_h, rhobar, c = bisection(alpha, rho0, m_h, rcut, rho_h, drhodr_h, mc_on_mh_min, mc_on_mh_max, eps_mc)
    m_c = m_h * c

    # Write out the final solution (uses a shorter dt, perhaps erroneously - need to check)
    xi, eta, theta = eval_rk(alpha, rho0, xi_h, rhobar, write=True)

    print '### Solution found ###'
    print 'rho0             = ', rho0
    print 'rho_h/theta**nnn = ', rho_h / theta[-1]**nnn
    print 'alpha            = ', alpha

    with open('out/param.out', 'w+') as open_file:
            open_file.write('{:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}\n'.format(p_dt, nnn, xi_h, xi[-1], rhobar, rho0, alpha, m_c))

    print "Completed run with no errors (we hope)!"


# Determines mass at the cut radius and at the lower boundary for bisection method
def mle_cut_mass(alpha, rho0, rhobar, rcut, rho_h, drhodr_h, m_h, a, c):

    # Perform iteration to find alpha and rho0 for specific value of m_c
    alpha, rho0, xi_h = newton_raphson(alpha, rho0, rcut, rhobar, rho_h, drhodr_h)

    # Get xi and theta to determine MLE mass
    xi, eta, theta = eval_rk(alpha, rho0, xi_h, rhobar)
    alphacgs = alpha * cm_Rsun
    mMLE = cumtrapz(4. * np.pi * alphacgs**2 * xi**2 * rho0 * theta**nnn, alphacgs * xi, initial=0.)

    m_c_a = a * m_h + (mMLE[-1]/g_Msun) # Mass at cut radius if we use lower limit of mc_on_mh
    m_c_c = c * m_h + (mMLE[-1]/g_Msun) # Mass at cut radius if we use middle value of mc_on_mh

    m_criterion_a = (m_c_a - m_h) / m_h # Percentage difference between cut mass and m_h
    m_criterion_c = (m_c_c - m_h) / m_h # if we use low and middle values of mc_on_mh, respectively

    return alpha, rho0, xi_h, m_criterion_a, m_criterion_c

# Bisection method for finding m_c such that, at the cut radius, mMLE + m_c = m_h
def bisection(alpha, rho0, m_h, rcut, rho_h, drhodr_h, a, b, tol):
    c = (a+b)/2.0 # Middle value of mc_on_mh
    while (b-a)/2.0 > tol:
        rhobar = 3. * m_h * c * g_Msun / (4. * np.pi * rcut**3 * cm_Rsun**3)

        print "\nMass of core particle is between {} Msun and {} Msun.".format(m_h * a, m_h * b)

        alpha, rho0, xi_h, m_criterion_a, m_criterion_c = mle_cut_mass(alpha, rho0, rhobar, rcut, rho_h, drhodr_h, m_h, a, c)

        if m_criterion_c == 0: # We're exactly right!
            return alpha, rho0, xi_h, rhobar, c
        elif m_criterion_a * m_criterion_c < 0: # If the low and middle criteria are of differing signs (i.e. they bracket the zero value)
            b = c                               # then bring the top mc_on_mh down to the middle value
        else:                                   # else
            a = c                               # bring the bottom value up to the top value
        c = (a+b)/2.0
    return alpha, rho0, xi_h, rhobar, c

if __name__ == '__main__':
    mle_run(kernel)
