# This script is made to be automatically ran every day to produce ocean maps in near real-time
# It downloads the latest data based on today's date from the CMEMS FTP server into the inputs folder,
# then runs the BFN-QG algorithm in MASSH to assimilate this data into a dynamical model and create SSH maps for every 3 hours.
# Finally, the program averages these values for each day and derives the geostrophic velocities and 
# normalized relative vorticity from the daily SSH, which it saves as .nc files.
# Lagragian diagnostics are run using the LOCEAN Lamta algorithm.
# Finally, everything is sent to external FTP servers for access by end users. 

import os
from datetime import timedelta
import sys


########################### Parameters to adjust ##########################################################################################

destination = None # Available options : 'ifremer',
make_lagrangian_diags = True # True or False

dir_massh = '/bettik/PROJECTS/pr-data-ocean/stellaa/MASSH/mapping'
path_config = './NRT_BFN_main_config.py' 


###########################################################################################################################################
###  0. INITIALIZATION
###########################################################################################################################################

sys.path.append(dir_massh)
currdir=os.getcwd()

from src import exp
config = exp.Exp(path_config)
today = config.EXP.final_date
numdays = int((today-config.EXP.init_date)/timedelta(days = 1))

lon_min = config.GRID.lon_min                            
lon_max = config.GRID.lon_max                               
lat_min = config.GRID.lat_min                                 
lat_max = config.GRID.lat_max
bbox = [lon_min, lon_max, lat_min, lat_max]   

from tools.plot_tools import where_is_this
where_is_this(bbox, 20)


###########################################################################################################################################
###  1. DATA DOWNLOAD
###########################################################################################################################################

from tools.ftp_transfer import download_nadirs_cmems, download_swot_nadir
from tools.remapping import make_mdt

# What datasets to download
datasets = [
    'dataset-duacs-nrt-global-al-phy-l3', 
    'dataset-duacs-nrt-global-c2n-phy-l3', 
    'dataset-duacs-nrt-global-h2b-phy-l3',
    'dataset-duacs-nrt-global-s3a-phy-l3',
    'dataset-duacs-nrt-global-s3b-phy-l3',
    'cmems_obs-sl_glo_phy-ssh_nrt_j3n-l3-duacs_PT1S',
    'cmems_obs-sl_glo_phy-ssh_nrt_s6a-hr-l3-duacs_PT1S',
]

dataset_l4 = 'dataset-duacs-nrt-global-merged-allsat-phy-l4'

# FTP connection to CMEMS server and data download
download_nadirs_cmems(currdir, today, numdays, datasets, dataset_l4)
download_swot_nadir(currdir, today)

# If needed, creates appropriate mdt file (mdt has to be downloaded already, though)
make_mdt(currdir,bbox)


############################################################################################################################################
### 2. BOUNDARY CONDITIONS
############################################################################################################################################

from tools.remapping import compute_filled_map

# Rework DUACS dataset for optimal boundary conditions : extrapolate data to fill coasts. 
# Then a mask is used in BFN to select only ocean and avoid awkward 0 values around coasts
BC_data_path = currdir+'/input/'+today.strftime('%Y%m%d')+'/dataset-duacs-nrt-global-merged-allsat-phy-l4/*.nc'
save_new_BC_to = currdir+'/input/'+today.strftime('%Y%m%d')+'/duacs_l4_filled.nc'

compute_filled_map(BC_data_path, save_new_BC_to, bbox)


############################################################################################################################################
### 3. DATA ASSIMILATION WITH BFN-QG
############################################################################################################################################

# State
from src import state as state
State = state.State(config)

# Obs
from src import obs as obs # if no files to open, re-download data
dict_obs = obs.Obs(config,State)

# Model
from src import mod as mod
Model = mod.Model(config,State)

# Bondary Conditions
from src import bc as bc
Bc = bc.Bc(config)

# Inversion
from src import inv as inv
inv.Inv(config,State,Model,dict_obs=dict_obs,Bc=Bc)


###########################################################################################################################################
### 4. RESULTS PROCESSING
###########################################################################################################################################

from tools.remapping import nc_processing
nc_processing(today=today)


#######################################################################################
### 5. LAMTA LAGRANGIAN DIAGNOSTICS
#######################################################################################

if make_lagrangian_diags == True:
    dir_lamta = '/bettik/PROJECTS/pr-data-ocean/stellaa/lamtaLR'
    from tools.remapping import apply_lamta
    lamta_diags_results = apply_lamta(currdir, dir_lamta, today, bbox, numdays=3, bathylvl =-1500)


###########################################################################################################################################
### 6. MAPS UPLOAD
###########################################################################################################################################
# Here, choose the right function to send to the right place. 

if destination == 'ifremer':
    from tools.ftp_transfer import ftp_to_ifremer
    ftp_to_ifremer(today, currdir)