import read_mesa
import write_bc
import mle_bc_nr
from parameters import *
import sys
import os

input_file = sys.argv[1]
output_filename = input_file.split('.')[0] + '.out'
param_filename = input_file.split('.')[0] + '.param'
output_file = 'out/{}'.format(output_filename)
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

alpha, xi_h, xi_max, rho0, rhobar, m_c, bc = mle_bc_nr.mle_run(input_file)

# Write out the final solution (uses a shorter dt, perhaps erroneously - need to check)
xi, eta, theta = mle_bc_nr.eval_rk(alpha, rho0, xi_h, rhobar, xi_max, write=True, filename=output_file)

rho_h = bc['rho']

print('\n### Solution found ###')
print('rho0             = {}'.format(rho0))
print('rho_h/theta**nnn = {}'.format(rho_h / theta[-1]**nnn))
print('alpha            = {}'.format(alpha))

print("Completed run with no errors (we hope)!")

with open('out/{}'.format(param_filename), 'w+') as open_file:
        open_file.write('{:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}  {:10.8E}\n'.format(p_dt, nnn, xi_h, xi[-1], rhobar, rho0, alpha, m_c))

