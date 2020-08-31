import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import sys

from .input_file_utils import read_input_file

# Reads a MESA data file, returns a SimpleNamespace containing the desired columns
def read_mesa_file(in_file):
    # Initialise list to store the data from the MESA file
    vals = []

    # Open MESA file
    with open(in_file, 'r') as open_file:

        # Loops over each line in the file
        for num, line in enumerate(open_file, 1):
            line_list = line.split()

            # Ignore the first 5 lines(which are the header)
            if num < 6:
                pass

            # Read line 6, which contains the headers
            elif num == 6:
                titles = line_list

            # Read all remaining lines into the vals array
            # Each column is stored in it's own list
            else:
                for i, val in enumerate(line_list):
                    try:
                        vals[i].append(float(val))
                    except IndexError:
                        vals.append([float(val)])

    # Specify which columns to retain from the MESA file
    cols_to_keep = [('mass', 'logM'), ('radius', 'logR'), ('pressure', 'logP'), ('rho', 'logRho'), ('energy', 'logE'), ('temperature', 'logT'), ('x_mass_fraction_H', ''), ('y_mass_fraction_He', ''), ('z_mass_fraction_metals', '')] 

    # Initialise dictionary to contain the index of the desired columns in the file
    output_cols = {}

    # Loop over the desired columns
    for cols in cols_to_keep:

        # Loop over all the column titles from the file
        for j, title in enumerate(titles):

            # If the column title is one of the desired columns
            # take note of the index of the column and whether or not it is logged
            if title in cols:
                col_logged = False
                if cols.index(title) == 1:
                    col_logged = True

                output_cols.update({cols[0] : (j,col_logged)})
                break
        
        # If the loop reaches the end, then the desired column is not in the file, and the run therefore cannot continue
        else:
            print("Neither {} nor {} found in the input profile!".format(*cols))
            if cols[0] == 'energy': # Energy is not strictly necessary, but it is interesting to see
                output_cols.update({cols[0]: None})
            else:
                exit()
    
    # Initialise dictionary to contain the column data
    output_dict = {}

    # Loop over the columns and their indices/logged state
    for col, col_data in output_cols.items():

        # If the col_data is None, the column was not in the file (at this point, only energy can do this)
        if col_data is None:
            output_vals = np.zeros(len(vals[0]))

        # Otherwise, use the column index to get the column data from the MESA file data
        else:
            output_vals = np.array(vals[col_data[0]])
        
        # If col_data[1] is True, the column was logged, which we don't want
        # Convert back
        if col_data is not None and col_data[1]:
            output_vals = 10**output_vals

        # Add the column data to the output dictionary, but flip it so that it increases in radius
        output_dict.update({col : np.flipud(output_vals)})

    # Aesthetic change
    # This is effectively the same as a dictionary, but I can get the items by,
    # e.g., output.radius rather than output['radius']
    output = SimpleNamespace(**output_dict)

    return output

# Reads a MESA data file including atomic species required by the EoS, returns a SimpleNamespace containing the desired columns
def read_mesa_file_detail(in_file):
    # Initialise list to store the data from the MESA file
    vals = []

    # Read in the input file from the command line
    profile_file = sys.argv[1]
    mle_input_file = profile_file.split('.')[0] + '.in'
    _, _, _, mesa_eos_params = read_input_file(mle_input_file)

    # Open MESA file
    with open(in_file, 'r') as open_file:

        # Loops over each line in the file
        for num, line in enumerate(open_file, 1):
            line_list = line.split()

            # Ignore the first 5 lines(which are the header)
            if num < 6:
                pass

            # Read line 6, which contains the headers
            elif num == 6:
                titles = line_list

            # Read all remaining lines into the vals array
            # Each column is stored in it's own list
            else:
                for i, val in enumerate(line_list):
                    try:
                        vals[i].append(float(val))
                    except IndexError:
                        vals.append([float(val)])

    # Specify which columns to retain from the MESA file
     #basic columns
    basic_cols_to_keep = [('mass', 'logM'), ('radius', 'logR'), ('pressure', 'logP'), ('rho', 'logRho'), ('energy', 'logE'), ('temperature', 'logT'), ('x_mass_fraction_H', ''), ('y_mass_fraction_He', ''), ('z_mass_fraction_metals', ''), ('abar', ''), ('zbar', '')] 
     #chemistry columns from input file, note that we have to remove the "i" character in front of the elements names in the input file
    chem_cols = [e[1:] for e in mesa_eos_params.chem_id]
    empty_cols = [''] * len(mesa_eos_params.chem_id)
    chem_cols_to_keep = list(zip(chem_cols,empty_cols))
     #merge lists
    cols_to_keep = basic_cols_to_keep + chem_cols_to_keep

    # Initialise dictionary to contain the index of the desired columns in the file
    output_cols = {}

    # Loop over the desired columns
    for cols in cols_to_keep:

        # Loop over all the column titles from the file
        for j, title in enumerate(titles):

            # If the column title is one of the desired columns
            # take note of the index of the column and whether or not it is logged
            if title in cols:
                col_logged = False
                if cols.index(title) == 1:
                    col_logged = True

                output_cols.update({cols[0] : (j,col_logged)})
                break
        
        # If the loop reaches the end, then the desired column is not in the file, and the run therefore cannot continue
        else:
            print("Neither {} nor {} found in the input profile!".format(*cols))
            if cols[0] == 'energy': # Energy is not strictly necessary, but it is interesting to see
                output_cols.update({cols[0]: None})
            else:
                exit()
    
    # Initialise dictionary to contain the column data
    output_dict = {}

    # Loop over the columns and their indices/logged state
    for col, col_data in output_cols.items():

        # If the col_data is None, the column was not in the file (at this point, only energy can do this)
        if col_data is None:
            output_vals = np.zeros(len(vals[0]))

        # Otherwise, use the column index to get the column data from the MESA file data
        else:
            output_vals = np.array(vals[col_data[0]])
        
        # If col_data[1] is True, the column was logged, which we don't want
        # Convert back
        if col_data is not None and col_data[1]:
            output_vals = 10**output_vals

        # Add the column data to the output dictionary, but flip it so that it increases in radius
        output_dict.update({col : np.flipud(output_vals)})

    # Aesthetic change
    # This is effectively the same as a dictionary, but I can get the items by,
    # e.g., output.radius rather than output['radius']
    output = SimpleNamespace(**output_dict)

    return output

# For reading the file supplied by Natasha Ivanova
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

    output_dict = {}

    for col, col_name in enumerate(output_cols):
        output_vals = np.array(vals[col])
        output_dict.update({output_cols[col] : output_vals})

    output = SimpleNamespace(**output_dict)

    return output

# Sets which function should be used to read the input data file
# This will typically be left on read_mesa_file
def read_mesa(in_file):
    return read_mesa_file(in_file)
def read_mesa_detail(in_file):
    return read_mesa_file_detail(in_file)

# Primarily used for testing functions
if __name__ == "__main__":
    x = read_ivanova_file("./ivanova_profile.data")
    y = read_mesa_file("./P12_MESA_Profile.data")

    print(x.radius,y.radius)
