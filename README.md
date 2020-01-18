# Modified Lane-Emden Equation Solver

This project is an implementation of the modified Lane-Emden equation method laid out by [Ohlmann et al. (2017)](https://www.aanda.org/articles/aa/abs/2017/03/aa29692-16/aa29692-16.html). It reconstructs the core of a stellar profile within a specified radius such that it is consistent with having the core condensed into a point mass with softened gravitational potential.

## Requirements
Python 3.7 or higher

pipenv (can be installed from the command line using `pip install pipenv`)

Git

## Installation

Clone a copy of this repo by running the command `git clone https://github.com/TomReichardt/ModifiedLE.git`

Change to the new ModifiedLE directory, `cd ModifiedLE`

Run `pipenv install` inside this directory, which should create a new virtual environment with all the dependencies from the Pipfile and Pipfile.lock

## Usage

Once all packages have been installed, we can open the virtual environment using `pipenv shell`

Copy your MESA profile into a new directory, e.g. ModifiedLE/profiles/ (or you can separate out your profiles even further with ModifiedLE/profiles/name_of_mesa_profile/)

Before we run the MLE solver, we need to create an input file. Luckily, there is a tool for that! From the base directory, run `python create_input_file.py profile`, where profile is the path to your MESA profile in the directory you created

`create_input_file.py` gives an interface to make initial guesses for the parameters alpha and rho0, along with being able to set the polytrope index n for the core region. Slide the sliders until the solution approximately lines up with the original profile at the cut radius.
