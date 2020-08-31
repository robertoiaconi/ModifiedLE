import numpy as np
from sys import float_info
from scipy.integrate import cumtrapz, solve_ivp
from types import SimpleNamespace
from copy import deepcopy
from .utils import phantom, ohlmann, g_Msun, cm_Rsun, cgs_G
from .bisection import bisection

# Initialise timestepping parameters
n_steps   = 1000    # number of points to evaluate solution at
tmin      = 1.e-6  # must start at finite value to avoid singularity

def modified_lane_emden(t, y, mle_params):
    theta, eta = y
    xi         = t

    u, chi, dchidu = mle_params.kernel.kernel(mle_params.alpha * xi)
    
    dfdt0 = -eta / xi**2
    dfdt1 = xi**2 * (theta**mle_params.n + mle_params.rhobar / mle_params.rho0 * (chi + 1./3. * u * dchidu))

    dfdt = np.array([dfdt0, dfdt1])

    return dfdt

def zero_crossing(t, y):
    return y[0]
zero_crossing.terminal = True
zero_crossing.direction = -1

def evaluate_mle(mle_params, write=False, filename='out/generic.out'):
    fun    = lambda t, y: modified_lane_emden(t, y, mle_params)
    t_span = (tmin, mle_params.xi_max)
    y0     = [1., 0.] #Initial conditions for integration
    t_eval = np.linspace(tmin, mle_params.xi_max, n_steps)

    sol = solve_ivp(fun=fun, t_span=t_span, y0=y0, t_eval=t_eval, events=zero_crossing, dense_output=True)

    xi    = np.array(sol.t)
    theta = np.array(sol.y[0])
    eta   = np.array(sol.y[1])

    if write:
        dthetadxi = -eta / xi**2
        chi, dchidu = np.array([]), np.array([])
        for ixi in xi:
            u, ichi, idchidu = mle_params.kernel.kernel(mle_params.alpha * ixi)
            chi = np.append(chi,ichi)
            dchidu = np.append(dchidu, idchidu)
        sol_out = np.array([xi, theta, dthetadxi, chi, dchidu])
        np.savetxt(filename, sol_out.T)

    return xi, eta, theta

def rhocut(mle_params, BCs):
    xi, eta, theta = evaluate_mle(mle_params)

    r = mle_params.alpha * xi
    ic = (np.abs(r-BCs.r)).argmin()
    rho = mle_params.rho0 * theta[ic]**mle_params.n
    rhop = mle_params.rho0 * theta[ic-1]**mle_params.n
    drhodr = (rho - rhop) / (r[ic] - r[ic-1])  #This could probably be improved
    return rho, drhodr

# Find how rho and drhodr vary as alpha varies around present value of alpha
def drhodalpha(mle_params, BCs):
    h_alpha = np.sqrt(float_info.epsilon) * mle_params.alpha
    rho, drho = [None]*4, [None]*4
    for i, val in enumerate([-2,-1,1,2]):
        mle_params_temp = deepcopy(mle_params)
        mle_params_temp.alpha += val * h_alpha
        rho[i], drho[i] = rhocut(mle_params_temp, BCs)

    d_rho = (rho[0] - 8.*rho[1] + 8.*rho[2] - rho[3]) / (12.*h_alpha)
    d_drho = (drho[0] - 8.*drho[1] + 8.*drho[2] - drho[3]) / (12.*h_alpha)
    return d_rho, d_drho

# Find how rho and drhodr vary as rho0 varies around present value of rho0
def drhodrho0(mle_params, BCs):
    h_rho0 = np.sqrt(float_info.epsilon) * mle_params.rho0
    rho, drho = [None]*4, [None]*4
    for i, val in enumerate([-2,-1,1,2]):
        mle_params_temp = deepcopy(mle_params)
        mle_params_temp.rho0 += val * h_rho0
        rho[i], drho[i] = rhocut(mle_params_temp, BCs)

    d_rho = (rho[0] - 8.*rho[1] + 8.*rho[2] - rho[3]) / (12.*h_rho0)
    d_drho = (drho[0] - 8.*drho[1] + 8.*drho[2] - drho[3]) / (12.*h_rho0)
    return d_rho, d_drho

# Newton-Raphson with two variables, matching rho = rho_h and drhodr = drhodr_h
def newton_raphson(mle_params, BCs):

    # Initialise Jacobian
    J = np.array([[0.,0.],[0.,0.]])

    for i in range(100): # 100 is arbitrary, but if it hasn't converged by then it likely won't

        mle_params.xi_h = mle_params.kernel.radius * mle_params.kernel.h / mle_params.alpha
        mle_params.xi_max = BCs.r / mle_params.alpha

        J[0][0], J[1][0] = drhodalpha(mle_params, BCs) # J = [drho/dalpha    d^2rho/drdalpha]
        J[0][1], J[1][1] = drhodrho0(mle_params, BCs)  #     [drho/drho0     d^2rho/drdrho0 ]

        J = np.linalg.inv(J) # Matrix inverse of J

        rho, drhodr = rhocut(mle_params, BCs)     # Find rho and drho/dr at the cut radius
        mle_params.alpha -= np.dot(J, [rho - BCs.rho, drhodr - BCs.drhodr])[0] # Increment alpha based on gradients
        mle_params.rho0 -= np.dot(J, [rho - BCs.rho, drhodr - BCs.drhodr])[1]   # Increment rho0 based on gradients

        rho_check = abs((BCs.rho - rho) / rho)
        drhodr_check = abs((BCs.drhodr - drhodr) / drhodr)

        print('\n##### Iteration {} #####'.format(i))
        print('alpha = {}'.format(mle_params.alpha))
        print('rho0  = {}'.format(mle_params.rho0))
        print('Desired cut values : rho = {}, drhodr = {}'.format(BCs.rho, BCs.drhodr))
        print('MLE cut values     : rho = {}, drhodr = {}'.format(rho, drhodr))
        print('Cut rho difference    = {}'.format(rho_check))
        print('Cut drhodr difference = {}'.format(drhodr_check))

        if (rho_check < 1e-9 and drhodr_check < 1e-9):
            break
    return mle_params


#Main function
def mle_run(mle_params, bisection_inputs, BCs):
    # Set softening length of the kernel
    mle_params.kernel.set_hsoft(BCs.r / mle_params.kernel.radius)

    # Guts of the code - use bisection root finder to iterate m_c until the cut mass matches m_h
    mc_on_mh, mle_params = bisection(bisection_inputs.mratio_min, bisection_inputs.mratio_max, bisection_inputs.tol,
                                     bisection_function, mle_params, BCs)

    # The mass of the core particle consistent with the calculate MLE profile
    mle_params.mc = BCs.m * mc_on_mh

    return mle_params


# Determines mass at the cut radius and at the lower boundary for bisection method
def mle_cut_mass(mle_params, BCs):
    # Perform iteration to find alpha and rho0 for specific value of m_c
    mle_params = newton_raphson(mle_params, BCs)

    # Get xi and theta to determine MLE mass
    xi, eta, theta = evaluate_mle(mle_params)
    alphacgs = mle_params.alpha * cm_Rsun
    mle_mass = cumtrapz(4. * np.pi * alphacgs**2 * xi**2 * mle_params.rho0 * theta**mle_params.n, alphacgs * xi, initial=0.) / g_Msun

    return mle_mass[-1], mle_params

def bisection_function(lower, upper, middle, mle_params, BCs):
    print("\nMass of core particle is between {} Msun and {} Msun.".format(BCs.m * lower, BCs.m * upper))
    mle_params.rhobar = 3. * BCs.m * middle * g_Msun / (4. * np.pi * BCs.r**3 * cm_Rsun**3)

    mle_mass, mle_params = mle_cut_mass(mle_params, BCs)

    lower_cut_mass = (lower * BCs.m) + mle_mass # Mass at cut radius if we use lower limit of mc_on_mh
    middle_cut_mass = (middle * BCs.m) + mle_mass # Mass at cut radius if we use middle value of mc_on_mh

    lower_mass_criterion = (lower_cut_mass - BCs.m) / BCs.m # Percentage difference between cut mass and m_h
    middle_mass_criterion = (middle_cut_mass - BCs.m) / BCs.m # if we use low and middle values of mc_on_mh, respectively

    return lower_mass_criterion, middle_mass_criterion, mle_params



if __name__ == '__main__':
    mle_run(kernel)
