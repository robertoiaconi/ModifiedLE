#python rewrite of read_mesa.pro
import numpy as np
import sys, os
from collections import namedtuple
import matplotlib.pyplot as plt

def read_mesa_file(in_file):
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

    cols_to_keep = [('mass', 'logM'), ('radius', 'logR'), ('pressure', 'logP'), ('rho', 'logRho'), ('energy', 'logE'), ('temperature', 'logT'), ('x_mass_fraction_H', ''), ('y_mass_fraction_He', ''), ('z_mass_fraction_metals', '')] 

    output_cols = {}

    for cols in cols_to_keep:
        for j, title in enumerate(titles):
            if title in cols:
                col_logged = False
                if cols.index(title) == 1:
                    col_logged = True

                output_cols.update({cols[0] : (j,col_logged)})
                break
        else:
            print("Neither {} nor {} found in the input profile!".format(*cols))
            if cols[0] == 'energy':
                output_cols.update({cols[0]: None})
            else:
                exit()
                
            #output_cols.update({cols[0] : (j,col_logged)})

    output_data = namedtuple('output_data', output_cols.keys())

    output_dict = {}

    for col, col_data in output_cols.items():
        if col_data is None:
            output_val = np.zeros(len(vals[0]))
        else:
            output_vals = np.array(vals[col_data[0]])
        
        if col_data is not None and col_data[1]:
            output_vals = 10**output_vals
        output_dict.update({col : output_vals})

    output = output_data(**output_dict)

    return output


def read_ivanova_file(in_file):
    vals = []

    with open(in_file, 'r') as open_file:
        for num, line in enumerate(open_file, 1):
            line_list = line.split()
            for i, ik in enumerate(line_list):
                try:
                    vals[i].append(float(ik))
                except IndexError:
                    vals.append([float(ik)])

    output_cols = ['radius', 'mass', 'entropy', 'energy', 'eth', 'enuc', 'enu', 'temperature', 'rho', 'pressure', 'ion_entropy', 'betab', 'x_mass_fraction_H', 'y_mass_fraction_He', 'z_mass_fraction_metals']

    output_data = namedtuple('output_data', output_cols)

    output_dict = {}

    for col, col_name in enumerate(output_cols):
        output_vals = np.flipud(np.array(vals[col]))
        output_dict.update({output_cols[col] : output_vals})

    output = output_data(**output_dict)

    return output

def read_mesa(in_file):
    return read_mesa_file(in_file)

if __name__ == "__main__":
    x = read_ivanova_file("./ivanova_profile.data")
    y = read_mesa_file("./P12_MESA_Profile.data")

    print(x.radius,y.radius)