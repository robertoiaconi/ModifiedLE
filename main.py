import read_mesa
import boundary_conditions
import modified_lane_emden
from parameters import *
import sys
import os
import json

input_file = sys.argv[1]
output_filename = input_file.split('.')[0] + '.out'
param_filename = input_file.split('.')[0] + '.param'
output_file = 'out/{}'.format(output_filename)
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

BCs = boundary_conditions.find_bc(input_file, cut_radius)
mle_params = modified_lane_emden.mle_run(kernel, BCs)

# Write out the final solution (uses a shorter dt, perhaps erroneously - need to check)
xi, eta, theta = modified_lane_emden.eval_rk(mle_params, write=True, filename=output_file)

print('\n### Solution found ###')
print('rho0             = {}'.format(mle_params.rho0))
print('rho_h/theta**nnn = {}'.format(BCs.rho / theta[-1]**mle_params.n))
print('alpha            = {}'.format(mle_params.alpha))

print('Core mass = {} Msun'.format(mle_params.mc))

print("Completed run with no errors (we hope)!")

with open('out/{}'.format(param_filename), 'w+') as open_file:
        open_file.write(json.dumps(vars(mle_params)))

