# Modified Lane-Emden Equation Solver

This project is an implementation of the modified Lane-Emden equation method laid out by [Ohlmann et al. (2017)](https://www.aanda.org/articles/aa/abs/2017/03/aa29692-16/aa29692-16.html). It reconstructs the core of a stellar profile within a specified radius such that it is consistent with having the core condensed into a point mass with softened gravitational potential.

## Requirements
Python 3.7 or higher

virtualenv

Git

## Installation

Clone a copy of this repo by running the command `git clone https://github.com/TomReichardt/ModifiedLE.git`

Change to the new ModifiedLE directory, `cd ModifiedLE`

Create a virtual environment with `virtualenv venv`

Start the new virtual environment with `venv\Scripts\activate.bat` (Windows) or `source venv/bin/activate` (macOS or Linux)

Once the virtual environment is running, use `pip install -r requirements.txt`

## Usage

Copy your MESA profile into a new directory within ModifiedLE/profiles/ (for example, I have a MESA profile called test_profile.data, hence I copy it into ModifiedLE/profiles/test_profile/)

### Creating an input file

We now need to produce an input file for the modified Lane-Emden solver. Use the command `python create_input_file.py profiles/test_profile/test_profile.data` (replacing `profiles/test_profile/test_profile.data` with the path to your profile)

This script opens a plot window with the input density profile and a modified Lane-Emden profile

Set rcut to the whatever value you want. The core will be reconstructed within this radius

Set the polytropic index, n, to whatever value (though typically 1.5 or 3 are used)

Slide the sliders for alpha and rho_0 to get the MLE profile to line up at the cut radius with the MESA profile

Press the 'Write file' button, and then close the plot window. You should have a .in file in the same directory as your MESA profile

You are now free to run the MLE solver!

### Running the MLE solver

You can now run the MLE solver with `python main.py profiles/test_profile/test_profile.data` (again, replace with the path to your profile)

With the .in file created in the previous section, this should quickly converge to the correct profile

### Creating the output file (and seeing the profile you have created)
