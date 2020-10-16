import sys
import json

from mle.boundary_conditions import find_bc
from mle.modified_lane_emden import mle_run, evaluate_mle
from mle.input_file_utils import write_input_file, read_input_file

profile_file = sys.argv[1]

output_file = profile_file.split('.')[0] + '.out'
input_file = profile_file.split('.')[0] + '.in'

mle_params, bisection_inputs, cut_radius, mesa_eos_params = read_input_file(input_file)

BCs = find_bc(profile_file, cut_radius)
mle_params = mle_run(mle_params, bisection_inputs, BCs)

# Obtain the final solution
xi, eta, theta = evaluate_mle(mle_params, write=True, filename=output_file)

print('\n### Solution found ###')
print('rho0             = {}'.format(mle_params.rho0))
print('rho_h/theta**nnn = {}'.format(BCs.rho / theta[-1]**mle_params.n))
print('alpha            = {}'.format(mle_params.alpha))

print('Core mass = {} Msun'.format(mle_params.mc))

print("Completed run with no errors (we hope)!")

# Write the final solution to a file
write_input_file(input_file, mle_params, bisection_inputs, cut_radius, mesa_eos_params)

