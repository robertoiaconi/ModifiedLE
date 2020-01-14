import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from modified_lane_emden import evaluate_mle
import boundary_conditions
from types import SimpleNamespace
from read_mesa import read_mesa
import utils
from utils import phantom, g_Msun, cm_Rsun, cgs_G

# This script is designed to help provide initial guesses for alpha and rho0 to the main program
# Slide the sliders until the gradient and value of the lines approximately match at the cut location
# Note the values of alpha and rho0, then edit parameters.py

def write_input_file(filename, mle_params, bisection_inputs, cut_radius):
    mle_params.kernel = str(mle_params.kernel)
    inputs = {
        'mle_params' : vars(mle_params),
        'bisection_inputs' : vars(bisection_inputs),
        'cut_radius' : cut_radius
    }
    with open(filename.split('.')[0] + '.in', 'w+') as open_file:
            json.dump(inputs, open_file, indent=4)

def read_input_file(filename):
    with open(filename) as open_file:
        inputs = json.load(open_file)

    mle_params = SimpleNamespace(**inputs['mle_params'])
    mle_params.kernel = getattr(utils, mle_params.kernel)
    bisection_inputs = SimpleNamespace(**inputs['bisection_inputs'])
    cut_radius = inputs['cut_radius']

    return mle_params, bisection_inputs, cut_radius

def main():
    # Bisection parameters
    bisection_inputs = SimpleNamespace(
        tol        = 0.00001, # sets tolerance on mc
        mratio_min = 0.2,    # lower guess for ratio of core mass to cut mass
        mratio_max = 0.99   # upper guess for ratio of core mass to cut mass
    )

    # Set the radius within which to soften the core, units of Rsun
    cut_radius = 1

    # The MESA profile you wish to core
    input_file = sys.argv[1]
    BCs = boundary_conditions.find_bc(input_file, cut_radius)

    mle_params = SimpleNamespace()

    # Initial parameters - values are chosen arbitrarily here to be approximately of the correct magnitude
    mle_params.alpha = 5.
    mle_params.rho0 = 0.005
    mle_params.n = 3

    # Set kernel
    mle_params.kernel = phantom
    mle_params.kernel.set_hsoft(BCs.r / mle_params.kernel.radius)

    # These initial parameters depend on the previous ones
    mle_params.xi_max = BCs.r / mle_params.alpha
    mle_params.mc = BCs.m
    mle_params.rhobar = 3. * mle_params.mc * g_Msun / (4. * np.pi * (mle_params.alpha * mle_params.xi_max)**3 * cm_Rsun**3)

    # Integrate the initial MLE equation
    xi, eta, theta = evaluate_mle(mle_params)

    # Convert xi and theta to radius and density 
    r = mle_params.alpha*xi
    rho = mle_params.rho0 * np.power(theta, mle_params.n)

    # Set up the plots
    fig, ax = plt.subplots()

    # Plot the input profile
    s = read_mesa(input_file)
    ax.loglog(s.radius, s.rho)

    # Plot the initial MLE solution
    l, = ax.loglog(r, rho)

    # Plot the cut radius as a vertical line
    ax.axvline(BCs.r, color="k", linestyle="-.")

    # Adjust the subplot to give room for sliders
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax.margins(x=0)

    # Create the slider axes
    axcolor = 'lightsteelblue'
    axalpha = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    axrho0 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    axn = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)

    # Create sliders
    salpha = Slider(axalpha, 'Alpha', 0.5, 20.0, valinit=mle_params.alpha, color='steelblue')
    srho0 = Slider(axrho0, 'Rho0', 0.0, 0.1, valinit=mle_params.rho0, color='steelblue')
    sn = Slider(axn, 'n', 0.0, 5.0, valinit=mle_params.n, valstep=0.5, color='steelblue')

    # Updater function for the sliders - maybe splitting this into three functions is better
    def update(val):
        mle_params.alpha = salpha.val
        mle_params.rho0 = srho0.val
        mle_params.n = sn.val
        mle_params.xi_max = BCs.r / mle_params.alpha
        mle_params.rhobar = 3. * mle_params.mc * g_Msun / (4. * np.pi * (mle_params.alpha * mle_params.xi_max)**3 * cm_Rsun**3)
        xi, eta, theta = evaluate_mle(mle_params)
        r = mle_params.alpha*xi
        rho = mle_params.rho0 * np.power(theta,mle_params.n)
        l.set_xdata(r)
        l.set_ydata(rho)
        fig.canvas.draw_idle()

    # Updating based on changes to the slider values
    salpha.on_changed(update)
    srho0.on_changed(update)
    sn.on_changed(update)

    # Write file button
    writeax = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(writeax, 'Write file', color=axcolor, hovercolor='0.975')
    click_func = lambda event: write_input_file(input_file, mle_params, bisection_inputs, cut_radius)
    button.on_clicked(click_func)

    # Self-explanatory
    plt.show()

if __name__ == "__main__":
    main()