import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, TextBox
from types import SimpleNamespace

from mle.modified_lane_emden import evaluate_mle
from mle.boundary_conditions import find_bc
from mle.read_mesa import read_mesa
from mle.utils import phantom, g_Msun, cm_Rsun, cgs_G
from mle.input_file_utils import write_input_file, read_input_file

# This script is designed to help provide initial guesses for alpha and rho0 to the main program
# Slide the sliders until the gradient and value of the lines approximately match at the cut location
# Note the values of alpha and rho0, then edit parameters.py

# Create state container
state = SimpleNamespace()

# Bisection parameters
state.bisection_inputs = SimpleNamespace(
    tol        = 0.00001, # sets tolerance on mc
    mratio_min = 0.2,    # lower guess for ratio of core mass to cut mass
    mratio_max = 0.99   # upper guess for ratio of core mass to cut mass
)

# Set the radius within which to soften the core, units of Rsun
state.cut_radius = 1

# The MESA profile you wish to core
input_file = sys.argv[1]
state.BCs = find_bc(input_file, state.cut_radius)

# Initial parameters - values are chosen arbitrarily here to be approximately of the correct magnitude
state.mle_params = SimpleNamespace()
state.mle_params.alpha = 5.
state.mle_params.rho0 = 0.005
state.mle_params.n = 3

#MESA EoS parameters - some default values to be updated by manually editing the input file
state.mesa_eos_params = SimpleNamespace()
state.mesa_eos_params.species = 22
state.mesa_eos_params.chem_id = ["ineut","ih1","iprot","ihe3","ihe4","ic12","in14","io16","ine20","img24","isi28","is32","iar36","ica40","iti44","icr48","icr60","ife52","ife54","ife56","ico56","ini56"]
state.mesa_eos_params.net_iso = ["ineut","ih1","iprot","ihe3","ihe4","ic12","in14","io16","ine20","img24","isi28","is32","iar36","ica40","iti44","icr48","icr60","ife52","ife54","ife56","ico56","ini56"]

# Set kernel
state.mle_params.kernel = phantom
state.mle_params.kernel.set_hsoft(state.BCs.r / state.mle_params.kernel.radius)

# These initial parameters depend on the previous ones
state.mle_params.xi_max = state.BCs.r / state.mle_params.alpha
state.mle_params.mc = state.BCs.m
state.mle_params.rhobar = 3. * state.mle_params.mc * g_Msun / (4. * np.pi * (state.mle_params.alpha * state.mle_params.xi_max)**3 * cm_Rsun**3)

# Integrate the initial MLE equation
xi, eta, theta = evaluate_mle(state.mle_params)

# Convert xi and theta to radius and density 
r = state.mle_params.alpha*xi
rho = state.mle_params.rho0 * np.power(theta, state.mle_params.n)

# Set up the plots
fig = plt.figure(figsize=(12,6), facecolor='lightgrey')
gs = fig.add_gridspec(ncols=2, nrows=20, width_ratios=[0.25,0.75])
plot_ax = fig.add_subplot(gs[:,1])

# Plot the input profile
s = read_mesa(input_file)
plot_ax.loglog(s.radius, s.rho, label='MESA Profile')

# Plot the initial MLE solution
l, = plot_ax.loglog(r, rho, label='MLE Solution')

# Plot the cut radius as a vertical line
line = plot_ax.axvline(state.BCs.r, color="k", linestyle="-.")

# Create the slider axes
axcolor = 'w'
axalpha = fig.add_subplot(gs[10,0], facecolor=axcolor, xscale='log')
axrho0 = fig.add_subplot(gs[12,0], facecolor=axcolor, xscale='log')
axn = fig.add_subplot(gs[8,0], facecolor=axcolor)

# Create sliders
salpha = Slider(axalpha, '$\\alpha$', 0.1, 100., valinit=state.mle_params.alpha, valfmt='%1.1e', color='steelblue')
srho0 = Slider(axrho0, '$\\rho_0$', 0.00001, 10.0, valinit=state.mle_params.rho0, valfmt='%1.1e', color='steelblue')
sn = Slider(axn, '$n$', 0.0, 5.0, valinit=state.mle_params.n, valstep=0.5, color='steelblue')

# Update plot based on new values
def update(mle_params):
    xi, _, theta = evaluate_mle(mle_params)
    r = mle_params.alpha*xi
    rho = mle_params.rho0 * np.power(theta,mle_params.n)
    l.set_xdata(r)
    l.set_ydata(rho)
    fig.canvas.draw_idle()

# Updater function for alpha
def update_alpha(val):
    state.mle_params.alpha = val
    state.mle_params.xi_max = state.BCs.r / state.mle_params.alpha
    state.mle_params.rhobar = 3. * state.mle_params.mc * g_Msun / (4. * np.pi * (state.mle_params.alpha * state.mle_params.xi_max)**3 * cm_Rsun**3)
    update(state.mle_params)

# Updater function for rho0
def update_rho0(val):
    state.mle_params.rho0 = val
    update(state.mle_params)

# Updater function for n
def update_n(val):
    state.mle_params.n = val
    update(state.mle_params)

# Updating based on changes to the slider values
salpha.on_changed(update_alpha)
srho0.on_changed(update_rho0)
sn.on_changed(update_n)

# Cut radius textbox
cutrax = fig.add_subplot(gs[6,0])
textbox = TextBox(cutrax, '$r_\mathrm{cut}$', initial=str(state.cut_radius), color=axcolor)

def update_cut_radius(val):
    try:
        state.cut_radius = float(val)
    except ValueError:
        print("Value in textbox must be a number")
        return
    state.BCs = find_bc(input_file, state.cut_radius)
    state.mle_params.kernel.set_hsoft(state.BCs.r / state.mle_params.kernel.radius)

    state.mle_params.xi_max = state.BCs.r / state.mle_params.alpha
    state.mle_params.mc = state.BCs.m
    state.mle_params.rhobar = 3. * state.mle_params.mc * g_Msun / (4. * np.pi * (state.mle_params.alpha * state.mle_params.xi_max)**3 * cm_Rsun**3)
    line.set_xdata(state.cut_radius)
    update(state.mle_params)

textbox.on_submit(update_cut_radius)

# Write file button
writeax = fig.add_subplot(gs[14,0])
button = Button(writeax, 'Write file', color=axcolor, hovercolor='0.975')
def click_func(event):
    state.mle_params.mc = None
    state.mle_params.xi_max = None
    state.mle_params.rhobar = None
    write_input_file(input_file, state.mle_params, state.bisection_inputs, state.cut_radius, state.mesa_eos_params)
button.on_clicked(click_func)

# Self-explanatory
handles, labels = plot_ax.get_legend_handles_labels()
plot_ax.legend(handles, labels, loc=1)
plt.show()
