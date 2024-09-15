# ISM property maps data release from TNG50

This repository contains all local ISM property maps generated from TNG50 simulations (https://www.tng-project.org/) for making predictions in post-processing using the Pressure-Regulated, Feedback-Modulated Star Formation Model (PRFM). For more details on the PRFM model (see Ostriker & Kim 2022  https://arxiv.org/abs/2206.00681). A full description of the equations, methods, and implementation guide can be found in Hassan et al. (2024).

# Repository Structure
  - **generate_ism_maps.py** Fast post-processing code using multiprocessing to generate various ISM property maps of average/midplane and projected densities within 1 kpc $\times$ 1 kpc columns of all selected galaxies from TNG50, matching the required resolution of TIGRESS simulations. The code can be used to generate these maps from any full snapshot of the TNG simulations (not limited to TNG50).
  - **utils.py** Includes utility functions used within the post-processing code.
  - **requirements.txt** Lists all dependencies required to run the code successfully. This file is generated using pipreqs (https://github.com/bndr/pipreqs).
  - **examples** Contains a Jupyter notebook (read_data_and_plot.ipynb) and some helper functions (user_functions.py) to demonstrate how the data can be read/used, along with clear steps to reproduce most of the plots in Hassan et al. (2024).

# Download data
The output of the post-processing code (**generate_ism_maps.py**) from TNG50 simulation at z=0 and z=2, stored as pickled pandas objects. [Click here to download](https://drive.google.com/drive/folders/1dsk6z7ugOJEedwfW6uhaDi8YYGvfPEF9?usp=sharing)



# Running the post-processing code

To run the code:

      $ python generate_ism_maps.py snap

snap refers to the snapshot ID, which can be found at https://www.tng-project.org/. Refer to the Snapshot Tables for details; for example, use snap=99 for z=0. The only required change in the code is to provide the full path to the TNG simulation directory, e.g., basePath='/Illustris_IllustrisTNG_public_data_release/L35n2160TNG/output/'. The code is flexible, allowing users to modify parameters such as pixel size and scales to compute midplane and projected quantities. These can be adjusted by changing the values for cell_size, int_scale, and mid_scale. Additionally, sample selection criteria are applied automatically but can be customized by modifying mask_sample to fit any desired conditions.

The code is highly efficient and utilizes multiprocessing to distribute computations across available cores. As a reference, the code processes 10,397 and 21,155 subhalos using 48 CPUs in 6.6 minutes and 14.6 minutes at z=0 and z=2, respectively.

# Output file structure
After downloading the data using the link above, each pickled pandas object (e.g., /data/all_data_z0.pkl or /data/all_data_z2.pkl) contains the following maps for each galaxy in the selected sample (see examples for how to read them)

| Column | Description |
| ----------- | ----------- |
| subhalo_id  | Subhalo ID from TNG50 subhalo catalogs |
| Sigma_g | Surface density of gas in $M_{\odot}/kpc^{2}$|
| rho_g | Volumetric density of gas in $M_{\odot}/kpc^{3}$|
| H_g | Measured scale height of gas in $kpc$ | 
|n_H | Hydrogen number density in $cm^{-3}$|
|Z_g| Mass weighted averagred metallicity in solar units|
|P_th| Thermal pressure in $K/cm^{3}$|
|P_mag|Magentic pressure in $K/cm^{3}$|
|P_turb| Turbulent pressure in $K/cm^{3}$ |
|Sigma_*| Surface density of stars in $M_{\odot}/kpc^{2}$ |
|rho_*|  Volumetric density of star in $M_{\odot}/kpc^{3}$|
|H_*| Stellar scale height in $kpc$|
|Sigma_SFR| Surface density of star formation in in $M_{\odot}/yr/kpc^{2}$|
|rho_dm| Volumetric density of dark matter in $M_{\odot}/kpc^{3}$|
