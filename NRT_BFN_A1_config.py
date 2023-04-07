#!/usr/bin/env python3
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 19:20:42 2021

@author: leguillou

modified by Adrien Stella
from config_2020a_BFNQG.py
"""

#################################################################################################################################
# Global libraries     
#################################################################################################################################
 
from datetime import datetime,timedelta,date
today = date.today().strftime('%Y%m%d')
#today = date.fromisoformat('2023-03-21').strftime('%Y%m%d')

numdays = 15
final_date = datetime.today() # or to stop at yday : - timedelta(days=1)
#final_date = date.fromisoformat('2023-03-21')
init_date = final_date - timedelta(days=numdays)
 
#################################################################################################################################
# EXPERIMENTAL PARAMETERS
#################################################################################################################################

name_experiment = 'NRT_BFN_A1'

EXP = dict(

    name_experiment = name_experiment, # name of the experiment

    saveoutputs = True, # save outputs flag (True or False)

    name_exp_save = name_experiment, # name of output files

    path_save = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/output/'+today+'/', # path of output files

    tmp_DA_path = "./scratch/NRT_BFN_A1", # temporary data assimilation directory path,

    init_date = init_date, # initial date (yyyy,mm,dd,hh) 

    final_date = final_date,  # final date (yyyy,mm,dd,hh) 

    assimilation_time_step = timedelta(hours=3),

    saveoutput_time_step = timedelta(hours=3),  # time step at which the states are saved 

    flag_plot = 0,

    write_obs = False, # save observation dictionary in *path_obs*

    path_obs = "./scratch/NRT_BFN_A1"

)
    
#################################################################################################################################
# GRID parameters
#################################################################################################################################
NAME_GRID = 'myGrid'

myGrid = dict(

    super = 'GRID_GEO',

    lon_min = -2.,                                        # domain min longitude

    lon_max = 10.,                                        # domain max longitude

    lat_min = 36.,                                         # domain min latitude

    lat_max = 44.,                                         # domain max latitude

    dlon = 1/16.,                                            # zonal grid spatial step (in degree)

    dlat = 1/16.,                                            # meridional grid spatial step (in degree)
 
    name_init_mask = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/duacs_l4_filled.nc',

    name_var_mask = {'lon':'longitude','lat':'latitude','var':'adt'} # adt variable is before filling. Different variable for BCs.

)

#################################################################################################################################
# Model parameters
#################################################################################################################################
NAME_MOD = 'myMOD'

myMOD = dict(

    super = 'MOD_QG1L_NP',

    name_var = {'SSH':"ssh", "PV":"pv"},

    dtmodel = 600, # model timestep

    c0 = 1.5, # phase speed of baroclinic 1st mode - calculé avec Emmanuel à partir de LR = 15km en Med.

)

#########################################
###   External boundary conditions    ### 
#########################################

NAME_BC = 'myBC' # For now, only BC_EXT is available

myBC = dict(

    super = 'BC_EXT',

    file = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/duacs_l4_filled.nc',

    name_lon = 'longitude',

    name_lat = 'latitude',

    name_time = 'time',

    name_var = {'SSH':'adt_full'}, # name of the boundary conditions variable

    name_mod_var = {'SSH':'ssh'},

    dist_sponge = 30 # Peripherical band width (km) on which the boundary conditions are applied

)

#################################################################################################################################
# Observation parameters
#################################################################################################################################
NAME_OBS = ['ALG','C2N','H2B','J3N','S3A','S3B','S6A']

ALG = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/dataset-duacs-nrt-europe-al-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

C2N = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/dataset-duacs-nrt-europe-c2n-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

H2B = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/dataset-duacs-nrt-europe-h2b-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

J3N = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/cmems_obs-sl_eur_phy-ssh_nrt_j3n-l3-duacs_PT0.2S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S3A = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/cmems_obs-sl_eur_phy-ssh_nrt_s3a-l3-duacs_PT0.2S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S3B = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/cmems_obs-sl_eur_phy-ssh_nrt_s3b-l3-duacs_PT0.2S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S6A = dict(

    super = 'OBS_SSH_NADIR',

    path = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today+'/cmems_obs-sl_eur_phy-ssh_nrt_s6a-hr-l3-duacs_PT0.2S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/cmems_obs-sl_med_phy-mdt_my_l4-0.0417deg_P20Y_1679318915395.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)


#################################################################################################################################
# INVERSION
#################################################################################################################################
NAME_INV = 'myINV'

myINV = dict(
    
    super = 'INV_BFN',

    window_size = timedelta(days=7), # length of the bfn time window

    window_output = timedelta(days=3), # length of the output time window, in the middle of the bfn window. (need to be smaller than *bfn_window_size*)

    propagation_timestep = timedelta(hours=3), # propagation time step of the BFN, corresponding to the time step at which the nudging term is computed

    criterion = 0.001, # convergence criterion. typical value: 0.01  

    max_iteration = 10, # maximal number of iterations if *bfn_criterion* is not met

)


#################################################################################################################################
# Diagnostics
#################################################################################################################################

NAME_DIAG = None