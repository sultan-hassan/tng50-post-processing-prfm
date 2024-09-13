# ISM property maps Data Release from TNG50

This repository contains all local ISM property maps generated from TNG50 simulations (https://www.tng-project.org/) to generate predictions in post-processing for the Pressure-Regulated, Feedback-Modulated Star Formation Model (PRFM), see Ostriker & Kim, https://arxiv.org/abs/2206.00681. For the full description of the equations/methods/implementation guide, please see Hassan+2024. 

# Repository Structure
  - **generate_ism_maps.py** Fast post-processing code using multiprocessing to generate various ISM property maps of average/midplane and projected densities within 1 kpc $\times$ 1 kpc columns of all selected galaxies from TNG50 to match the required resolution of TIGRESS simulations. The code can be used to generate all these maps from any full snapshot of the TNG simualtions (not specifically for TNG50).  
  - **utils.py** includes some utility functions that are used within the prost-processing code.  
  - **requirements.txt** includes all dependencies required to run the code sucessfully. This file is generated using pipreqs (https://github.com/bndr/pipreqs).
  - **data** contains the output of the post-processing code from TNG50 at z=0 and z-2, which are stored as pickled pandas objects.
  - **examples** contains a jupyter notebook (read_data_and_plot.ipynb) and some helper functions (user_functions.py) to show examples of the data can be used and provide clear steps to reproduce most of the plots in Hassan+2024.

# Running the post-processing code

To tun the code:

  python generate_ism_maps.py snap

snap is the snapshot id, which comes from https://www.tng-project.org/, see Snap tables, e.g. use snap=99 for z=0.

