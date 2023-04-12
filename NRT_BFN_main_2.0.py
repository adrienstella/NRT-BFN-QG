# This script is made to be automatically ran every day at 6am as an OAR job by submit_NRT_BFN.ksh.
# It downloads the latest data based on today's date from the CMEMS FTP server, into the inputs folder,
# then runs the BFN-QG algorithm in MASSH to assimilate this data and output ssh maps for every 3 hours.
# Finally, the program averages these values for each day and derives the geostrophic velocities and 
# normalized relative vorticity from the daily ssh, which it saves as .nc files that are also sent to 
# Ifremer's FTP server. 

from ftplib import FTP
import numpy as np
import os
import pandas as pd
from datetime import datetime, timedelta, date
import xarray as xr
import pyinterp.backends.xarray
import pyinterp.fill
import netCDF4

###########################################################################################################################################
###  1. DOWNLOAD DATA
###########################################################################################################################################

today = date.today()
# today = date.fromisoformat('2023-03-21') # To get data from another day
numdays = 15

# Define spatial domain (!) should be better coded to use it from config.
lon_min = -2                              
lon_max = 10                                
lat_min = 36                                  
lat_max = 44                                    

# What data to download and where to put it
datasets = [
    'dataset-duacs-nrt-europe-al-phy-l3', 
    'dataset-duacs-nrt-europe-c2n-phy-l3', 
    'dataset-duacs-nrt-europe-h2b-phy-l3',
    'cmems_obs-sl_eur_phy-ssh_nrt_j3n-l3-duacs_PT0.2S',
    'cmems_obs-sl_eur_phy-ssh_nrt_s3a-l3-duacs_PT0.2S',
    'cmems_obs-sl_eur_phy-ssh_nrt_s3b-l3-duacs_PT0.2S',
    'cmems_obs-sl_eur_phy-ssh_nrt_s6a-hr-l3-duacs_PT0.2S',
]

inputs_location='/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'

# FTP connection to CMEMS server and data download
os.chdir('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/')
currdir=os.getcwd()

# Set user name and password
username = 'username'
password = 'password'

# Connect to the ftp server
ftp = FTP('nrt.cmems-du.eu',username,password)

# Choose dates to download
first_day = today - timedelta(days=numdays)
dates = pd.date_range(end = today, periods = numdays).to_pydatetime().tolist()

from tools.ftp_transfer import ftp_cmems_download_month

# Download all L3 products
for i in np.arange(0,len(datasets)):

    os.makedirs(inputs_location+today.strftime('%Y%m%d')+'/'+datasets[i], exist_ok = True)
    os.chdir(inputs_location+today.strftime('%Y%m%d')+'/'+datasets[i])

    ftp_cmems_download_month(ftp, '/Core/SEALEVEL_EUR_PHY_L3_NRT_OBSERVATIONS_008_059/', datasets[i], str(today.year), str(today.month))
    os.chdir(inputs_location+today.strftime('%Y%m%d')+'/'+datasets[i])

    if (today.month != first_day.month): # (!) only works when the assimilation duration doesn't span over more than two months
            ftp_cmems_download_month(ftp, '/Core/SEALEVEL_EUR_PHY_L3_NRT_OBSERVATIONS_008_059/', datasets[i], str(first_day.year), str(first_day.month))
            os.chdir(inputs_location+today.strftime('%Y%m%d')+'/'+datasets[i])

os.chdir(currdir)
print('Obs data downloaded successfully')


# Download DUACS L4 product
dataset_l4 = 'dataset-duacs-nrt-europe-merged-allsat-phy-l4'
print('Retreiving data for dataset '+dataset_l4)

os.makedirs(inputs_location+today.strftime('%Y%m%d')+'/'+dataset_l4, exist_ok = True)
os.chdir(inputs_location+today.strftime('%Y%m%d')+'/'+dataset_l4)

ftp_cmems_download_month(ftp, '/Core/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/', dataset_l4, str(today.year), str(today.month))
os.chdir(inputs_location+today.strftime('%Y%m%d')+'/'+dataset_l4)

if (today.month != first_day.month): # (!) only works when the assimilation duration doesn't span over more than two months
    ftp_cmems_download_month(ftp, '/Core/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/', dataset_l4, str(first_day.year), str(first_day.month))

os.chdir(currdir)
print('DUACS L4 data downloaded successfully')

ftp.quit()

############################################################################################################################################
### 2. COMPUTE BOUNDARY CONDITIONS
############################################################################################################################################

# Rework the DUACS dataset for optimal boundary conditions : extrapolate data to fill coasts. 
# Then a mask is used in BFN to select only ocean and avoid awkward 0 values around coasts

ds = xr.open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today.strftime('%Y%m%d')+'/dataset-duacs-nrt-europe-merged-allsat-phy-l4/*.nc')
ds = ds.sel(longitude = slice(lon_min,lon_max), latitude = slice(lat_min,lat_max))

longitude = ds.longitude.values
latitude = ds.latitude.values

x_axis = pyinterp.Axis(longitude)
y_axis = pyinterp.Axis(latitude)

adt_filled = np.zeros((len(ds.time), len(longitude), len(latitude)))
adt_filled_upright = np.zeros((len(ds.time), len(latitude), len(longitude)))

for t in np.arange(len(ds.time)):
    has_converged, adt_filled[t] = pyinterp.fill.gauss_seidel(pyinterp.Grid2D(x_axis, y_axis, ds.adt[t].values.T)) # values are transposed you have to T them back
    adt_filled_upright[t] = adt_filled[t].T

ds['adt_full'] = (['time', 'latitude', 'longitude'], adt_filled_upright)

ds.to_netcdf('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/input/'+today.strftime('%Y%m%d')+'/duacs_l4_filled.nc', mode = 'w')


############################################################################################################################################
### 3. RUN DATA ASSIMILATION WITH BFN-QG
############################################################################################################################################

# Config
path_config = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/NRT_BFN_A1_config.py'  

dir_massh = '/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/MASSH-A1/mapping'
import sys
sys.path.append(dir_massh)

from src import exp
config = exp.Exp(path_config)

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
### 4. PROCESS RESULTS
###########################################################################################################################################

# Load the last week and take daily averages
bfn_output = xr.open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/output/'+today.strftime('%Y%m%d')+'/*.nc', concat_dim='time', combine='nested')
bfn_output_dates = bfn_output.assign(time=pd.to_datetime(bfn_output.time.dt.date))
last_week = bfn_output_dates.where(bfn_output_dates.time >= pd.to_datetime(today-timedelta(days=6)), drop =True)
daily_mean_ssh = last_week.groupby("time").mean("time")

# A few physical parameters for computations
g=9.81 # m/s2
c0=1.5, # phase speed of baroclinic 1st mode - calculé avec Emmanuel à partir de LR = 15km en Med.
earth_rad=6.371e6 # 6371km
earth_w=2*np.pi/86400 # rad/s

# Create the dy, dx and f grids based on physical parameters
lon=daily_mean_ssh.lon.values
lat=daily_mean_ssh.lat.values
x, y = np.meshgrid(lon, lat)

dy = np.ones(y.shape)
dx = np.ones(x.shape)

dy[1:-1,:]=earth_rad*2*np.pi/360*(y[1:-1,:]-y[0:-2,:])
dx[:,1:-1]=earth_rad*2*np.pi/360*(x[:,1:-1]-x[:,0:-2])*np.cos(((y[:,1:-1]+y[:,0:-2])/2)*np.pi/180)

# Neumann condition at boundaries
dy[0,:]=dy[1,:]
dy[-1,:]=dy[-2,:]
dx[:,0]=dx[:,1]
dx[:,-1]=dx[:,-2]

f=2*earth_w*np.sin(y*np.pi/180)


# U, V & PV computation
u=np.zeros((len(daily_mean_ssh.time), len(lat), len(lon)))
v=np.zeros((len(daily_mean_ssh.time), len(lat), len(lon)))
xi_norm=np.zeros((len(daily_mean_ssh.time), len(lat), len(lon)))

from tools.vars import h2uv, h2rv

for t in np.arange(len(daily_mean_ssh.time)):
    (u[t, :, :], v[t, :, :]) = h2uv(ssh = daily_mean_ssh.ssh[t, :, :], dy = dy, dx = dx, g = g, f = f)
    xi_norm[t, :, :] = h2rv(ssh = daily_mean_ssh.ssh[t, :, :], dy = dy, dx = dx, g = g, f = f)

daily_mean_ssh['u'] = (['time', 'lat', 'lon'],  u)
daily_mean_ssh['v'] = (['time', 'lat', 'lon'],  v)
daily_mean_ssh['xi_norm'] = (['time', 'lat', 'lon'],  xi_norm)


# Save final netcdf files 
os.makedirs('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/maps/'+today.strftime('%Y%m%d')+'/', exist_ok = True)
for d in daily_mean_ssh.time.values:
    date = pd.to_datetime(d)
    ds=daily_mean_ssh.where(daily_mean_ssh.time == d, drop=True)
    ds.to_netcdf('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/maps/'+today.strftime('%Y%m%d')+'/NRT_BFN_'+date.strftime('%Y%m%d')+'.nc', mode = 'w')




###########################################################################################################################################
### 5. UPLOAD MAPS
###########################################################################################################################################

os.chdir('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/')
currdir=os.getcwd()

# Set user name and password
username = 'username'
password = 'password'

# Connect to the ftp server
ftp = FTP('ftp.ifremer.fr',username,password)
ftp.cwd('MEDSSH_BFN')

ftp.mkd(today.strftime('%Y%m%d'))
ftp.cwd(today.strftime('%Y%m%d'))

# put the file in the directory on ftp
for d in pd.to_datetime(daily_mean_ssh.time.values):
    ftp.storbinary('STOR BFN_CSWOT_'+d.strftime('%Y%m%d')+'.nc', open('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT_BFN/maps/'+today.strftime('%Y%m%d')+'/NRT_BFN_'+d.strftime('%Y%m%d')+'.nc', 'rb'))

print('Files in the directory :')
ftp.nlst()
print('Closing connexion :')
ftp.quit()