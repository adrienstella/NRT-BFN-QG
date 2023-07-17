# NRT-BFN-QG

A wrapper for F. Le Guillou's BFN-QG algorithm which produces dynamic maps of ocean surface topography by assimilation of altimetric data into a quasi-geostrophic model. The present program presents a simplified, fully automatized workflow, usable with a crontab, hence permitting near-real-time applications. 

In particular, the following features are integrated here: data download, pre-processing, visualization of input data, output processing, validation of output data, lagrangian diagnostics, maps upload to external FTP server. When run at the regional scale, this algorithm produces high resolution SSH and geostrophic currents maps, resolving much finer dynamical structures than DUACS, which can be used for better informed glider deployment and ship routing, for example.

See example/ for an example notebook with figures.

### Main steps of the algorithm: 
1) Data download:
- Download alongtrack (L3) altimetric data from a constellation of 7 nadirs on the CMEMS FTP server
- Dowload SWOT nadir from AVISO
- Download gridded (L4) DUACS ADT product on the CMEMS FTP server (for boundary conditions & mask)
- Download MDT on the CMEMS FTP server

2) Pre-processing of MDT and boundary conditions for optimal use. 

3) Running the Back-and-Forth Nudging algorithm from MASSH (https://github.com/leguillf/MASSH) with the specified configuration: assimilating the data with the quasi-geostrophic model with the chosen parameters.
- Optional: plot input data

4) Outputs processing : cut off spin-up period, make daily averages, compute additional variables, save new .nc files.

5) Optional: Plot validation curves by comparing final maps and DUACS L4 maps to observation tracks.

6) Optional: Make lagrangian diagnostics with LAMTA and save associated results and plots.

7) Optional: Copy final .nc files to remote server (e.g., Ifremer FTP)

### Default final files:
- Daily NetCDF files in folders named after 'final_date', in separate folders named after config files' 'name_experiment'.
- 7 files for each run to match the 7 days of the latest assimilation window (can be adjusted for reanalysis mode).
- Variables inclued : Sea Surface Height, Geostrophic (horizontal) velocities, Normalized relative vorticity, Potential vorticity.

### Main parameters to set: 

##### In the configuration file:
- *final_date*: final date for assimilation (default = current computer day and time, for NRT. Can be changed for reanalysis mode)
- *numdays*: number of days of assimilation to run before final date
- *lon_min*, *lon_max*, *lat_min*, *lat_max*: spatial domain
- *dlon*, *dlat*: spatial resolution (recommended 1/16°, or 1/8° for a large domain)
- *c0*: phase speed of baroclinic first mode (can be found here for example: https://ceoas.oregonstate.edu/rossby_radius)

##### In the main script/notebook:
- *destination*: where to send the output maps (e.g. FTP)
- *make_lagrangian_diags*: whether or not to perform LAMTA analysis (requires access to private repo https://github.com/rousseletL/lamtaLR)
- *draw_L3*: whether or not to plot the input observation data tracks
- *make_alongtrack_rmse*: whether or not to plot the validation of the final maps against input data and DUACS

- *dir_massh*: location of MASSH/mapping folder
- *path_config*: location of the configuration file to use for this experiment

### Installation
:computer: _**How to get started ?**_

(a) Clone this repo: 
```
git clone https://github.com/adrienstella/NRT-BFN-QG.git
```
(b) Clone the MASSH repo in the same place:
```
git clone https://github.com/leguillf/MASSH.git
```
(c) Check paths match your configuration in the main script

(d) Enter your usernames and passwords for the FTP servers (create your own secretcodes.py file with your cmems / aviso accounts).

(e) Create your own experiment by making a copy of the NRT_BFN_main_config.py file and adjusting the parameters to fit your needs. Don't forget to point to that new config file in the main script's 'path_config' variable.

You're good to go! Just run NRT_BFN_main.py (or the notebook version) to start mapping.