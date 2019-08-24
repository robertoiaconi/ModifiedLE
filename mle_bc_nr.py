import os
import sys
import numpy as np
from scipy.integrate import cumtrapz, solve_ivp
from utils import phantom, ohlmann
import write_bc
from parameters import *
import matplotlib.pyplot as plt


def modified_lane_emden(t, y, alpha, xi_h, nnn, rhobar, rho0):
    theta, eta = y
    xi         = t

    dthetadxi = -eta/(xi**2)

    u, chi, dchidu = kernel.kernel(alpha * xi)
    
    dfdt0 = -eta / xi**2
    dfdt1 = xi**2 * (theta**nnn + rhobar / rho0 * (chi + 1./3. * u * dchidu))

    dfdt = np.array([dfdt0, dfdt1])

    return dfdt

def zero_crossing(t, y):
    return y[0]
zero_crossing.terminal = True
zero_crossing.direction = -1

def eval_rk(alpha, rho0, xi_h, rhobar, xi_max, write=False, filename='out/generic.out'):
    fun    = lambda t, y: modified_lane_emden(t, y, alpha, xi_h, nnn, rhobar, rho0)
    t_span = (tmin, xi_max)
    y0     = [1., 0.] #Initial conditions for integration
    t_eval = np.linspace(tmin, xi_max, 1000)

    sol = solve_ivp(fun=fun, t_span=t_span, y0=y0, t_eval=t_eval, events=zero_crossing, dense_output=True)

    xi    = np.array(sol.t)
    theta = np.array(sol.y[0])
    eta   = np.array(sol.y[1])

    if write:
        dthetadxi = -eta / xi**2
        chi, dchidu = np.array([]), np.array([])
        for ixi in xi:
            u, ichi, idchidu = kernel.kernel(alpha * ixi)
            chi = np.append(chi,ichi)
            dchidu = np.append(dchidu, idchidu)
        sol_out = np.array([xi, theta, dthetadxi, chi, dchidu])
        np.savetxt(filename, sol_out.T)

    return xi, eta, theta

def rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha):
    xi, eta, theta = eval_rk(alpha, rho0, xi_h, rhobar, xi_max)

    r = alpha * xi
    ircut = (np.abs(r-rcut)).argmin()
    rho = rho0 * theta[ircut]**nnn
    rhop = rho0 * theta[ircut-1]**nnn
    drhodr = (rho - rhop) / (r[ircut] - r[ircut-1])  #This could probably be improved
    return rho, drhodr

# Find how rho and drhodr vary as alpha varies around present value of alpha
def drhodalpha(rcut, rhobar, nnn, rho0, alpha, xi_h, xi_max):
    h_alpha = np.sqrt(sys.float_info.epsilon) * alpha
    rhon2, drhon2 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha-2*h_alpha)
    rhon1, drhon1 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha-h_alpha)
    rhop1, drhop1 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha+h_alpha)
    rhop2, drhop2 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha+2*h_alpha)
    drho = (rhon2 - 8.*rhon1 + 8.*rhop1 - rhop2) / (12.*h_alpha)
    ddrho = (drhon2 - 8.*drhon1 + 8.*drhop1 - drhop2) / (12.*h_alpha)
    return drho, ddrho

# Find how rho and drhodr vary as rho0 varies around present value of rho0
def drhodrho0(rcut, rhobar, nnn, rho0, alpha, xi_h, xi_max):
    h_rho0 = np.sqrt(sys.float_info.epsilon) * rho0
    rhon2, drhon2 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0-2*h_rho0, alpha)
    rhon1, drhon1 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0-h_rho0,   alpha)
    rhop1, drhop1 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0+h_rho0,   alpha)
    rhop2, drhop2 = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0+2*h_rho0, alpha)
    drho = (rhon2 - 8.*rhon1 + 8.*rhop1 - rhop2) / (12.*h_rho0)
    ddrho = (drhon2 - 8.*drhon1 + 8.*drhop1 - drhop2) / (12.*h_rho0)
    return drho, ddrho

# Newton-Raphson with two variables, matching rho = rho_h and drhodr = drhodr_h
def newton_raphson(alpha, rho0, rcut, rhobar, rho_h, drhodr_h):

    # Initialise Jacobian
    J = np.array([[0.,0.],[0.,0.]])

    for i in range(100): # 100 is arbitrary, but if it hasn't converged by then it likely won't

        xi_h = kernel.radius * kernel.h / alpha
        xi_max = rcut / alpha

        print(xi_h, xi_max)
        J[0][0], J[1][0] = drhodalpha(rcut, rhobar, nnn, rho0, alpha, xi_h, xi_max) # J = [drho/dalpha    d^2rho/drdalpha]
        J[0][1], J[1][1] = drhodrho0(rcut, rhobar, nnn, rho0, alpha, xi_h, xi_max)  #     [drho/drho0     d^2rho/drdrho0 ]

        J = np.linalg.inv(J) # Matrix inverse of J

        rho, drhodr = rhocut(rcut, rhobar, nnn, xi_h, xi_max, rho0, alpha)     # Find rho and drho/dr at the cut radius
        alpha = alpha - np.dot(J, [rho - rho_h, drhodr - drhodr_h])[0] # Increment alpha based on gradients
        rho0  = rho0  - np.dot(J, [rho - rho_h, drhodr - drhodr_h])[1] # Increment rho0 based on gradients

        rho_check = abs((rho_h-rho) / rho)
        drhodr_check = abs((drhodr_h - drhodr) / drhodr)

        print('\n##### Iteration {} #####'.format(i))
        print('alpha = {}'.format(alpha))
        print('rho0  = {}'.format(rho0))
        print('Desired cut values : rho = {}, drhodr = {}'.format(rho_h, drhodr_h))
        print('MLE cut values     : rho = {}, drhodr = {}'.format(rho, drhodr))
        print('Cut rho difference    = {}'.format(rho_check))
        print('Cut drhodr difference = {}'.format(drhodr_check))

        if (rho_check < 1e-9 and drhodr_check < 1e-9):
            break
    return alpha, rho0, xi_h, xi_max


#Main function
def mle_run(input_file):
    # Initialise boundary conditions and initial guesses for alpha and rho0
    alpha = alpha_init
    rho0  = rho0_init

    # Set softening length of the kernel
    kernel.set_hsoft(input_h)

    # Guts of the code - use bisection root finder to iterate m_c until the cut mass matches m_h
    alpha, rho0, xi_h, xi_max, rhobar, c, bc = bisection(input_file, alpha, rho0, min_mh, max_mh, eps_mc)

    return alpha, xi_h, xi_max, rho0, rhobar, c, bc


# Determines mass at the cut radius and at the lower boundary for bisection method
def mle_cut_mass(alpha, rho0, rhobar, bc, a, c):

    # Perform iteration to find alpha and rho0 for specific value of m_c
    alpha, rho0, xi_h, xi_max = newton_raphson(alpha, rho0, bc['r'], rhobar, bc['rho'], bc['drhodr'])

    # Get xi and theta to determine MLE mass
    xi, eta, theta = eval_rk(alpha, rho0, xi_h, rhobar, xi_max)
    alphacgs = alpha * cm_Rsun
    mMLE = cumtrapz(4. * np.pi * alphacgs**2 * xi**2 * rho0 * theta**nnn, alphacgs * xi, initial=0.)

    m_at_cut = input_mc + (mMLE[-1]/g_Msun) # Mass at cut radius if we use lower limit of mc_on_mh
                                                                                                   #CHANGE THESE LINES USING A AND C
    m_criterion_a = (a - m_at_cut) / m_at_cut # Percentage difference between cut mass and m_h
    m_criterion_c = (c - m_at_cut) / m_at_cut # if we use low and middle values of mc_on_mh, respectively

    return alpha, rho0, xi_h, xi_max, m_criterion_a, m_criterion_c

# Bisection method for finding m_c such that, at the cut radius, mMLE + m_c = m_h
def bisection(input_file, alpha, rho0, a, b, tol):
    c = (a+b)/2.0 # Middle value of mc_on_mh
    while (b-a)/2.0 > tol:
        bc = write_bc.find_bc(input_file, c)
        print(bc)
        rhobar = 3. * input_mc * g_Msun / (4. * np.pi * bc['r']**3 * cm_Rsun**3)

        print("a = {}, b = {}, c = {}".format(a, b, c))

        alpha, rho0, xi_h, xi_max, m_criterion_a, m_criterion_c = mle_cut_mass(alpha, rho0, rhobar, bc, a, c)

        if m_criterion_c == 0: # We're exactly right!
            return alpha, rho0, xi_h, rhobar, c
        elif m_criterion_a * m_criterion_c < 0: # If the low and middle criteria are of differing signs (i.e. they bracket the zero value)
            b = c                               # then bring the top mc_on_mh down to the middle value
        else:                                   # else
            a = c                               # bring the bottom value up to the top value
        c = (a+b)/2.0
    return alpha, rho0, xi_h, xi_max, rhobar, c, bc

if __name__ == '__main__':
    mle_run(kernel)
