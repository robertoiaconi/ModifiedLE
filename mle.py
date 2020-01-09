import os
import sys

import numpy as np
from scipy.integrate import cumtrapz, solve_ivp

import write_bc
from lesolver import LaneEmdenSolver
from utils import phantom, ohlmann
from parameters import *

# Main function
def mle_run(kernel, BCs):
    # Initialise solver
    solver = LaneEmdenSolver(alpha_init, rho0_init, nnn)
    
    # Set the kernel for the solver
    kernel.set_hsoft(BCs.r / kernel.radius)
    solver.set_kernel(kernel)

    # Set the boundary conditions
    solver.set_boundary_conditions(BCs)

    # Guts of the code - use bisection root finder to iterate m_c until the cut mass matches m_h
    bisection(solver, mc_on_mh_min, mc_on_mh_max, eps_mc)
    
    m_c = BCs.m * c

    return alpha, xi_h, xi_max, rho0, rhobar, m_c

# Bisection method for finding m_c such that, at the cut radius, mMLE + m_c = m_h
def bisection(solver, lower, upper, tol):
    middle = 0.5 * (upper + lower) # Middle value of mc_on_mh

    while 0.5 * (upper - lower) > tol:
        # Sets the mass of the core particle along with rhobar
        solver.set_core_mass(middle)

        print("\nMass of core particle is between {} Msun and {} Msun.".format(BCs.m * lower / g_Msun, BCs.m * upper / g_Msun))


        ### Need to fix below here ###
        alpha, rho0, xi_h, xi_max, m_criterion_a, m_criterion_c = mle_cut_mass(alpha, rho0, rhobar, BCs, a, c)

        if m_criterion_c == 0: # We're exactly right!
            return alpha, rho0, xi_h, rhobar, c
        elif m_criterion_a * m_criterion_c < 0: # If the low and middle criteria are of differing signs (i.e. they bracket the zero value)
            b = c                               # then bring the top mc_on_mh down to the middle value
        else:                                   # else
            a = c                               # bring the bottom value up to the top value
        c = (a+b)/2.0
    return alpha, rho0, xi_h, xi_max, rhobar, c


# Newton-Raphson with two variables, matching rho = rho_h and drhodr = drhodr_h
def newton_raphson(alpha, rho0, rhobar, BCs):

    # Initialise Jacobian
    J = np.array([[0.,0.],[0.,0.]])

    for i in range(100): # 100 is arbitrary, but if it hasn't converged by then it likely won't

        xi_h = kernel.radius * kernel.h / alpha
        xi_max = BCs.r / alpha

        J[0][0], J[1][0] = drhodalpha(BCs.r, rhobar, nnn, rho0, alpha, xi_h, xi_max) # J = [drho/dalpha    d^2rho/drdalpha]
        J[0][1], J[1][1] = drhodrho0(BCs.r, rhobar, nnn, rho0, alpha, xi_h, xi_max)  #     [drho/drho0     d^2rho/drdrho0 ]

        J = np.linalg.inv(J) # Matrix inverse of J

        rho, drhodr = rhocut(BCs.r, rhobar, nnn, xi_h, xi_max, rho0, alpha)     # Find rho and drho/dr at the cut radius
        alpha = alpha - np.dot(J, [rho - BCs.rho, drhodr - BCs.drhodr])[0] # Increment alpha based on gradients
        rho0 = rho0 - np.dot(J, [rho - BCs.rho, drhodr - BCs.drhodr])[1]   # Increment rho0 based on gradients

        rho_check = abs((BCs.rho-rho) / rho)
        drhodr_check = abs((BCs.drhodr - drhodr) / drhodr)

        print('\n##### Iteration {} #####'.format(i))
        print('alpha = {}'.format(alpha))
        print('rho0  = {}'.format(rho0))
        print('Desired cut values : rho = {}, drhodr = {}'.format(BCs.rho, BCs.drhodr))
        print('MLE cut values     : rho = {}, drhodr = {}'.format(rho, drhodr))
        print('Cut rho difference    = {}'.format(rho_check))
        print('Cut drhodr difference = {}'.format(drhodr_check))

        if (rho_check < 1e-9 and drhodr_check < 1e-9):
            break
    return alpha, rho0, xi_h, xi_max


# Determines mass at the cut radius and at the lower boundary for bisection method
def mle_cut_mass(solver, a, c):

    # Perform iteration to find alpha and rho0 for specific value of m_c
    # FIX THIS alpha, rho0, xi_h, xi_max = newton_raphson(alpha, rho0, rhobar, BCs)

    # Get xi and theta to determine MLE mass
    solver.update_profile()
    alphacgs = alpha * cm_Rsun
    mMLE = cumtrapz(4. * np.pi * alphacgs**2 * xi**2 * rho0 * theta**nnn, alphacgs * xi, initial=0.)

    m_c_a = a * BCs.m + (mMLE[-1]/g_Msun) # Mass at cut radius if we use lower limit of mc_on_mh
    m_c_c = c * BCs.m + (mMLE[-1]/g_Msun) # Mass at cut radius if we use middle value of mc_on_mh

    m_criterion_a = (m_c_a - BCs.m) / BCs.m # Percentage difference between cut mass and m_h
    m_criterion_c = (m_c_c - BCs.m) / BCs.m # if we use low and middle values of mc_on_mh, respectively

    return alpha, rho0, xi_h, xi_max, m_criterion_a, m_criterion_c



if __name__ == '__main__':
    mle_run(kernel)
