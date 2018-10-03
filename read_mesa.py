#python rewrite of read_mesa.pro
import numpy as np
import sys, os
from collections import namedtuple
import matplotlib.pyplot as plt

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

    cols_to_keep = [('mass', 'logM'), ('radius', 'logR'), ('pressure', 'logP'), ('rho', 'logRho'), ('energy', 'logE')] 

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
            print "Neither {} nor {} found in the input profile!".format(*cols)
            exit()
            #output_cols.update({cols[0] : (j,col_logged)})

    output_data = namedtuple('output_data', output_cols.keys())

    output_dict = {}

    for col, col_data in output_cols.items():
        output_vals = np.array(vals[col_data[0]])
        if col_data[1]:
            output_vals = 10**output_vals
        output_dict.update({col : output_vals})

    output = output_data(**output_dict)

    return output

