# Modified Lane-Emden Equation Solver

This project is an implementation of the modified Lane-Emden equation method laid out by [Ohlmann et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017A%26A...599A...5O/abstract). It reconstructs the core of a stellar profile within a specified radius such that it is consistent with having the core condensed into a point mass with softened gravitational potential. This version of the code uses the tabulated equation of state from the 1D stellar evolution code MESA by [Paxton et a. (2011)](https://ui.adsabs.harvard.edu/abs/2011ApJS..192....3P/abstract) to reconstruct the necessary thermodynamic quantities.

## Requirements
Python 3.7 or higher

virtualenv

Git

gfort2py (https://github.com/rjfarmer/gfort2py)

MESA

pyMesa (https://github.com/rjfarmer/pyMesa)

## Installation

Clone a copy of this repo by running the command `git clone https://github.com/robertoiaconi/ModifiedLE_MESAEoS.git`

Change to the new ModifiedLE directory, `cd ModifiedLE_MESAEoS`

Create a virtual environment with `virtualenv venv`

Start the new virtual environment with `venv\Scripts\activate.bat` (Windows) or `source venv/bin/activate` (macOS or Linux)

Once the virtual environment is running, use `pip install -r requirements.txt`

Install MESA, gfort2py and pyMesa, in this order, by follwing the instructions written in the respective websites

## Usage

Copy your MESA profile into a new directory within ModifiedLE_MESAEoS/profiles/ (for example, I have a MESA profile called test_profile.data, hence I copy it into ModifiedLE_MESAEoS/profiles/test_profile/)

### Creating an input file

We now need to produce an input file for the modified Lane-Emden solver. Use the command `python create_input_file.py profiles/test_profile/test_profile.data` (replacing `profiles/test_profile/test_profile.data` with the path to your profile)

This script opens a plot window with the input density profile and a modified Lane-Emden profile

![alt text][input1]

Set rcut to the whatever value you want. The core will be reconstructed within this radius

Set the polytropic index, n, to whatever value (though typically 1.5 or 3 are used)

Slide the sliders for &alpha; and &rho;<sub>0</sub> to get the MLE profile to line up at the cut radius with the MESA profile

![alt text][input2]

Press the 'Write file' button, and then close the plot window. You should have a .in file in the same directory as your MESA profile, which, in this case, looks like this:
```
{
    "mle_params": {
        "alpha": 4.363650352681207,
        "rho0": 0.015575837407081461,
        "n": 3.0,
        "kernel": "phantom",
        "xi_max": null,
        "mc": null,
        "rhobar": null,
    },
    "bisection_inputs": {
        "tol": 1e-05,
        "mratio_min": 0.2,
        "mratio_max": 0.99
    },
    "cut_radius": 8.0,
    "mesa_eos_params": {
        "species": 22,
        "chem_id": [
            "ineut",
            "ih1",
            "iprot",
            "ihe3",
            "ihe4",
            "ic12",
            "in14",
            "io16",
            "ine20",
            "img24",
            "isi28",
            "is32",
            "iar36",
            "ica40",
            "iti44",
            "icr48",
            "icr60",
            "ife52",
            "ife54",
            "ife56",
            "ico56",
            "ini56"
        ],
        "net_iso": [
            "ineut",
            "ih1",
            "iprot",
            "ihe3",
            "ihe4",
            "ic12",
            "in14",
            "io16",
            "ine20",
            "img24",
            "isi28",
            "is32",
            "iar36",
            "ica40",
            "iti44",
            "icr48",
            "icr60",
            "ife52",
            "ife54",
            "ife56",
            "ico56",
            "ini56"
        ]
    }
}
```

The `null`s are intentionally there to avoid potential confusion, as they are calculated within the main script

All the MESA equation of state parameters, i.e., `species`, `chem_id`s and `net_iso`s, are default values and you should replace them with those in your profile, in this case the file test_profile.data, before moving to the next step

You are now free to run the MLE solver!

### Running the MLE solver

You can now run the MLE solver with `python main.py profiles/test_profile/test_profile.data` (again, replace with the path to your profile)

With the .in file created in the previous section, this should quickly converge to the correct profile, with the final output looking something like

```
### Solution found ###
rho0             = 0.07235131076621393
rho_h/theta**nnn = 0.07235131076621397
alpha            = 1.2335651838892818
Core mass = 0.2577200441208392 Msun
Completed run with no errors (we hope)!
```

A .out file will have been created in the same profile directory, which contains the MLE profile variables (xi, theta, eta, etc.) while the .in file will have had the `null`s replaced with the calculated values

### Creating the output file (and seeing the profile you have created)

Finally, we can create the output profile with the command 
```python final_profile.py profiles/test_profile/test_profile.data profiles/test_profile/test_profile_cored.data```
where the second path will be the location and name of the reconstructed MLE profile

Upon running this script, a window should appear, showing you a myriad of values for the newly created profile

![alt text][output]

Admittedly, there could be some work done to improve this screen, but hopefully you get the idea

The script will have produced a new profile at the specified location, which looks something like

```
# This profile has been constructed to have a point mass core of m_c = 7.269 Msun, and a softening length of h_soft = 4.0 Rsun.
# [    Mass    ]  [  Pressure  ]  [Temperature ]  [   Radius   ]  [  Density   ]  [   E_int    ]
  1.97944099E+34  4.27375421E+02  3.38829688E+03  5.56718957E+13  1.97218305E-09  1.82378231E+12
  1.97944099E+34  4.27375422E+02  3.38829688E+03  5.56718957E+13  1.97218305E-09  1.82378231E+12
  1.97944099E+34  4.27375424E+02  3.38829688E+03  5.56718957E+13  1.97218307E-09  1.82378231E+12
  1.97944099E+34  4.27375433E+02  3.38829688E+03  5.56718957E+13  1.97218310E-09  1.82378231E+12
  1.97944099E+34  4.27375460E+02  3.38829688E+03  5.56718957E+13  1.97218321E-09  1.82378231E+12
  1.97944099E+34  4.27375554E+02  3.38829688E+03  5.56718956E+13  1.97218359E-09  1.82378231E+12
  1.97944099E+34  4.27375747E+02  3.38829710E+03  5.56718955E+13  1.97218437E-09  1.82378233E+12
  1.97944099E+34  4.27376087E+02  3.38829740E+03  5.56718954E+13  1.97218574E-09  1.82378236E+12
  1.97944099E+34  4.27376638E+02  3.38829796E+03  5.56718951E+13  1.97218796E-09  1.82378242E+12
  1.97944098E+34  4.27377697E+02  3.38829902E+03  5.56718946E+13  1.97219223E-09  1.82378254E+12
  ...
```
where all units are CGS except for the point mass values at the top

Happy simulating!

[input1]: https://github.com/TomReichardt/ModifiedLE/blob/simplified/assets/input_tool_1.png "Adjust parameters as you see fit"
[input2]: https://github.com/TomReichardt/ModifiedLE/blob/simplified/assets/input_tool_2.png "Looks pretty close"
[output]: https://github.com/TomReichardt/ModifiedLE/blob/simplified/assets/output_tool_1.png "Look, it's not my best work"
