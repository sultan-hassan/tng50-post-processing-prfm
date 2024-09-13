# ISM property maps Data Release from TNG50

This repository contains all local ISM property maps generated from TNG50 simulations (https://www.tng-project.org/) to generate predictions in post-processing for the Pressure-Regulated, Feedback-Modulated Star Formation Model (PRFM), see Ostriker & Kim, https://arxiv.org/abs/2206.00681. For the full description of the equations/methods/implementation guide, please see Hassan+2024. 

# Repository Structure
  - **generate_ism_maps.py** Fast post-processing code using multiprocessing to generate various ISM property maps of average/midplane and projected densities within 1 kpc $\times$ 1 kpc columns of all selected galaxies from TNG50 to match the required resolution of TIGRESS simulations. The code can be used to generate all these maps from any full snapshot of the TNG simualtions (not specifically for TNG50).  
  - **utils.py** includes some utility functions that are used within the prost-processing code.  
  - **requirements.txt** includes all dependencies required to run the code sucessfully. This file is generated using pipreqs (https://github.com/bndr/pipreqs).
  - **data** contains the output of the post-processing code from TNG50 at z=0 and z-2, which are stored as pickled pandas objects.
  - **examples** contains a jupyter notebook (read_data_and_plot.ipynb) and some helper functions (user_functions.py) to show examples of the data can be used and provide clear steps to reproduce most of the plots in Hassan+2024.

# Running the post-processing code

To run the code:

      $ python generate_ism_maps.py snap

snap is the snapshot id, which comes from https://www.tng-project.org/, see Snap tables, e.g. use snap=99 for z=0. The only change inside the code that is require is to provide the full path for TNG simulation directory e.g. basePath='/Illustris_IllustrisTNG_public_data_release/L35n2160TNG/output/'. The code is flexible to allow the user to change any parameters to any values to generate maps with different pixel size, and different scales to compute the midplane and projected quantities (e.g. simply change values for cell_size, int_scale, mid_scale). The code automatically applies the sample selection criteria, which can also be changed by modifying *mask_sample* to any conditions. 

The code is very efficient and uses multiprocessing to distribute computing to different available cores. As a reference, the code processes 10397 and 21155 subhalos using 48 cpus in 6.6 mins and 14.6 mins at z=0 and z=2, respectively.

# Output file structure (data directory)

| Quantity    | Description |
| ----------- | ----------- |
| subhalo_id  | Subhalo ID from TNG50 subhalo catalogs |
| Sigma_g | Surface density of gas in units of $M_{\odot}/kpc^{3}$|
