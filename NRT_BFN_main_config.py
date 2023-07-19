#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 19:20:42 2021

@author: leguillou

modified by Adrien Stella on 2023-07-19 from the MASSH example config_2020a_BFNQG.py
"""

#################################################################################################################################
# General Parameters    
#################################################################################################################################
 
from datetime import datetime,timedelta,date,time

final_date = datetime.combine(date.today(), time())
#final_date = datetime.combine(date.fromisoformat('2023-06-13'), time())
numdays = 8

today = final_date.strftime('%Y%m%d')
init_date = final_date - timedelta(days=numdays)
 
#################################################################################################################################
# EXPERIMENTAL PARAMETERS
#################################################################################################################################

name_experiment = 'NRT_BFN_main_config' 

EXP = dict(

    name_experiment = name_experiment, # name of the experiment

    saveoutputs = True, # save outputs flag (True or False)

    name_exp_save = name_experiment, # name of output files

    path_save = './output_'+name_experiment+'/'+today+'/', # path of output files

    tmp_DA_path = "./scratch/"+name_experiment, # temporary data assimilation directory path,

    init_date = init_date, # initial date (yyyy,mm,dd,hh) 

    final_date = final_date,  # final date (yyyy,mm,dd,hh) 

    assimilation_time_step = timedelta(hours=3),

    saveoutput_time_step = timedelta(hours=3),  # time step at which the states are saved 

    flag_plot = 1,

    write_obs = False, # save observation dictionary in *path_obs*

    path_obs = "./scratch/"+name_experiment

)

#################################################################################################################################
# GRID parameters
#################################################################################################################################
NAME_GRID = 'myGrid'

myGrid = dict(

    super = 'GRID_GEO',

    lon_min = -62.,                                        # domain min longitude

    lon_max = -52.,                                        # domain max longitude

    lat_min = 6.,                                         # domain min latitude

    lat_max = 17.,                                         # domain max latitude

    dlon = 1/8.,                                            # zonal grid spatial step (in degree)

    dlat = 1/8.,                                            # meridional grid spatial step (in degree)
 
    name_init_mask = './input_'+name_experiment+'/'+today+'/duacs_l4_filled.nc',

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

    c0 = 2.2, # phase speed of baroclinic 1st mode (m/s)

    dist_sponge_bc = 30

)

#########################################
###   External boundary conditions    ### 
#########################################

NAME_BC = 'myBC' # For now, only BC_EXT is available

myBC = dict(

    super = 'BC_EXT',

    file = './input_'+name_experiment+'/'+today+'/duacs_l4_filled.nc',

    name_lon = 'longitude',

    name_lat = 'latitude',

    name_time = 'time',

    name_var = {'SSH':'adt_full'}, # name of the boundary conditions variable

    name_mod_var = {'SSH':'ssh'}

)

#################################################################################################################################
# Observation parameters
#################################################################################################################################
NAME_OBS = ['ALG','C2N','H2B','J3N','S3A','S3B','S6A','SWOTN']

ALG = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/dataset-duacs-nrt-global-al-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

C2N = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/dataset-duacs-nrt-global-c2n-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

H2B = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/dataset-duacs-nrt-global-h2b-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

J3N = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/cmems_obs-sl_glo_phy-ssh_nrt_j3n-l3-duacs_PT1S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S3A = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/dataset-duacs-nrt-global-s3a-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S3B = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/dataset-duacs-nrt-global-s3b-phy-l3/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

S6A = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/cmems_obs-sl_glo_phy-ssh_nrt_s6a-hr-l3-duacs_PT1S/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(days=1)},

)

SWOTN = dict(

    super = 'OBS_SSH_NADIR',

    path = './input_'+name_experiment+'/'+today+'/nrt_global_swonc_phy_l3_1hz/',

    name_time = 'time',
    
    name_lon = 'longitude',

    name_lat = 'latitude',
    
    name_var = {'SSH':'sla_unfiltered'},

    add_mdt = True,

    path_mdt = './input_'+name_experiment+'/cnes_mdt_local.nc',

    name_var_mdt = {'lon':'longitude','lat':'latitude','mdt':'mdt'},
    
    nudging_params_ssh = {'sigma':0,'K':0.7,'Tau':timedelta(hours=10)},

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