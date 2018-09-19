#python rewrite of read_mesa.pro
import numpy as np
import sys, os
from collections import namedtuple

def read_mesa(in_file):
    vals = []

    with open(in_file, 'r') as open_file:
        for num, line in enumerate(open_file, 1):
            line_list = line.split()
            if num < 6:
                pass
            elif num == 6:
                titles = line_list
            else:
                for i, ik in enumerate(line_list):
                    try:
                        vals[i].append(float(ik))
                    except IndexError:
                        vals.append([float(ik)])

    cols_to_keep = ['mass', 'radius', 'pressure', 'pgas_div_ptotal', 'logRho', 'energy'] 

    cols = []

    for col in cols_to_keep:
        for j, title in enumerate(titles):
            if col == title:
                cols.append(j)

    output_data = namedtuple('output_data', cols_to_keep)

    output_vals = []

    for i, icol in enumerate(cols):
        output_vals.append(vals[icol])

    output_vals = np.array(output_vals)

    output_data.mass = output_vals[0]
    output_data.radius = output_vals[1]
    output_data.pressure = output_vals[2]
    output_data.pgas_div_ptotal = output_vals[3]
    output_data.logRho = output_vals[4]
    output_data.energy = output_vals[5]

    return output_data

