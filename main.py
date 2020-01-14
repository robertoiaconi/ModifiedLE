import read_mesa
import boundary_conditions
import modified_lane_emden
from input_file_utils import write_input_file, read_input_file
import sys
import os
import json

profile_file = sys.argv[1]

output_file = profile_file.split('.')[0] + '.out'
input_file = profile_file.split('.')[0] + '.in'

mle_params, bisection_inputs, cut_radius = read_input_file(input_file)

BCs = boundary_conditions.find_bc(profile_file, cut_radius)
mle_params = modified_lane_emden.mle_run(mle_params, bisection_inputs, BCs)

# Write out the final solution (uses a shorter dt, perhaps erroneously - need to check)
xi, eta, theta = modified_lane_emden.evaluate_mle(mle_params, write=True, filename=output_file)

print('\n### Solution found ###')
print('rho0             = {}'.format(mle_params.rho0))
print('rho_h/theta**nnn = {}'.format(BCs.rho / theta[-1]**mle_params.n))
print('alpha            = {}'.format(mle_params.alpha))

print('Core mass = {} Msun'.format(mle_params.mc))

print("Completed run with no errors (we hope)!")

write_input_file(input_file, mle_params, bisection_inputs, cut_radius)

