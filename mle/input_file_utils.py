import json
from types import SimpleNamespace
from mle import utils

def write_input_file(filename, mle_params, bisection_inputs, cut_radius, mesa_eos_params):
    mle_params.kernel = str(mle_params.kernel)
    inputs = {
        'mle_params' : vars(mle_params),
        'bisection_inputs' : vars(bisection_inputs),
        'cut_radius' : cut_radius,
	'mesa_eos_params' : vars(mesa_eos_params)
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
    mesa_eos_params = SimpleNamespace(**inputs['mesa_eos_params'])

    return mle_params, bisection_inputs, cut_radius, mesa_eos_params
