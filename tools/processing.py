import xarray as xr
import numpy as np
import pyinterp.backends.xarray
import pyinterp.fill
from pandas import to_datetime
from datetime import timedelta
import os

# Check MDT is ready
def make_mdt(name_experiment, currdir, bbox, dataset_mdt = "mdt_hybrid_cnes_cls18_cmems2020_global.nc"):

    """
    Checks if the Mean Dynamic Topography (MDT) file is available for a given experiment, and creates it if needed.

    Args:
        name_experiment (str): The name of the experiment.
        currdir (str): The current directory.
        bbox (list): The bounding box coordinates [lon_min, lon_max, lat_min, lat_max].
        dataset_mdt (str): The dataset MDT filename.

    Returns:
        None

    The function checks if the local MDT file exists. If not, it downloads the MDT file and creates a local copy
    within the experiment directory. The MDT is then subsetted to the specified bounding box.
    """

    if(not(os.path.isfile(currdir+'/input_'+name_experiment+'/cnes_mdt_local.nc'))):
        from tools.ftp_transfer import download_mdt
        download_mdt(name_experiment, currdir, dataset_mdt)
        ds = xr.open_dataset(currdir+'/input_'+name_experiment+'/'+dataset_mdt)
        ds = ds.sel(longitude = slice(bbox[0],bbox[1]), latitude = slice(bbox[2],bbox[3]))
        mdt = ds.mdt
        mdt.to_netcdf(currdir+'/input_'+name_experiment+'/cnes_mdt_local.nc')

def compute_filled_map(BC_data_path, save_to, bbox):

    """
    Computes and saves the filled map based on the provided boundary condition data.

    Args:
        BC_data_path (str): The path to the boundary condition data.
        save_to (str): The path to save the filled map.
        bbox (list): The bounding box coordinates [lon_min, lon_max, lat_min, lat_max].

    Returns:
        None

    The function opens the boundary condition data and subsets it to the specified bounding box. It then computes the
    filled map using the Gaussian Seidel method and saves the result to the specified location. The point of this operation
    is that the filled map can be masked by a different resolution land mask without having ill-valued pixels at coastlines.
    """

    ds = xr.open_mfdataset(BC_data_path)
    ds = ds.sel(longitude = slice(bbox[0],bbox[1]), latitude = slice(bbox[2],bbox[3]))

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

    ds.to_netcdf(save_to, mode = 'w')


def compute_u_v_rv(ssh_map):
        
    """
    Computes U, V, and RV (Relative Vorticity) based on the given sea surface height map.

    Args:
        ssh_map (xarray.Dataset): Sea surface height map containing the 'ssh' variable with dimensions (time, lat, lon).

    Returns:
        xarray.Dataset: The input dataset with additional variables 'u', 'v', and 'xi_norm' computed based on the SSH map.

    The function computes the U (eastward velocity), V (northward velocity), and RV (Relative Vorticity) fields using the
    sea surface height (SSH) map. It calculates the necessary physical parameters, such as dy, dx, and f, based on the
    geographical grid of the SSH map. Then, it iterates over each time step in the SSH map to compute the U, V, and RV fields.

    The computed U, V, and RV fields are added as additional variables 'u', 'v', and 'xi_norm' to the input SSH map dataset,
    respectively. The resulting dataset is returned.

    Note:
        - The 'ssh_map' dataset must have the 'ssh' variable with dimensions (time, lat, lon).
    """

    # A few physical parameters for computations
    g=9.81 # m/s2
    earth_rad=6.371e6 # 6371km
    earth_w=2*np.pi/86400 # rad/s

    # Create the dy, dx and f grids based on physical parameters
    lon=ssh_map.lon.values
    lat=ssh_map.lat.values
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

    # U, V & RV computation
    u=np.zeros((len(ssh_map.time), len(lat), len(lon)))
    v=np.zeros((len(ssh_map.time), len(lat), len(lon)))
    xi_norm=np.zeros((len(ssh_map.time), len(lat), len(lon)))

    from tools.vars import h2uv, h2rv

    for t in np.arange(len(ssh_map.time)):
        (u[t, :, :], v[t, :, :]) = h2uv(ssh = ssh_map.ssh[t, :, :], dy = dy, dx = dx, g = g, f = f)
        xi_norm[t, :, :] = h2rv(ssh = ssh_map.ssh[t, :, :], dy = dy, dx = dx, g = g, f = f)

    ssh_map['u'] = (['time', 'lat', 'lon'],  u)
    ssh_map['v'] = (['time', 'lat', 'lon'],  v)
    ssh_map['xi_norm'] = (['time', 'lat', 'lon'],  xi_norm)

    return ssh_map


def nc_processing(name_experiment, today, numdays = 6, frequency_hours = 24):

    """
    Performs netCDF processing for the given experiment and time period.

    Args:
        name_experiment (str): The name of the experiment.
        today (datetime.date): The current date.
        numdays (int, optional): The number of days to consider for processing. Default is 6.
        frequency_hours (int, optional): The frequency in hours for time averaging. Default is 24.

    Returns:
        None

    The function performs netCDF processing for the specified experiment and time period. It opens the output files
    from the experiment, calculates the mean sea surface height (SSH), and computes additional variables such as U, V,
    and PV. The processed data is saved as netCDF files.

    If `frequency_hours` is set to 24 (default), the function calculates the daily mean SSH by grouping the data
    by day and taking the mean. If `frequency_hours` is set to a different value, the function manually defines time
    bins of the specified frequency and calculates the mean SSH within each bin. The resulting time bins are labeled
    with the midpoints of each bin.
    """
    

    # Load the interest period and make binned averages
    bfn_output = xr.open_mfdataset('./output_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/*.nc', concat_dim='time', combine='nested')
    interest_period = bfn_output.where(bfn_output.time >= to_datetime(today-timedelta(days=numdays)), drop =True)

    interest_period = compute_u_v_rv(interest_period)

    if frequency_hours == 24:
        interest_period = interest_period.assign(time=to_datetime(interest_period.time.dt.date))
        binned_ssh = interest_period.groupby("time").mean("time")
    else:
        from pandas import date_range
        #Manually define the bins we want to use for time averaging
        time_averaging_bins = date_range(start=to_datetime(today-timedelta(days=numdays)), end=today, freq=timedelta(hours=frequency_hours))
        
        #Group our DataArray by these bins
        time_midpoints = date_range(start=to_datetime(today-timedelta(days=numdays))+timedelta(hours=frequency_hours)/2, end=today, freq=timedelta(hours=frequency_hours))
        binned_ssh = interest_period.groupby_bins(group='time', bins=time_averaging_bins, labels=time_midpoints).mean('time')
        binned_ssh = binned_ssh.rename({'time_bins': 'time'})


    # Save final netcdf files - this can be adapted depending on the needs of the receiver (variable names, how many files, etc.)
    os.makedirs('./maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/', exist_ok = True)
    os.makedirs('./maps_'+name_experiment+'/full_timeseries/', exist_ok = True)

    if frequency_hours == 24:
        format = '%Y%m%d'
    else:
        format = '%Y%m%d-%H'
        
    for d in binned_ssh.time.values:
        date = to_datetime(d)
        ds=binned_ssh.where(binned_ssh.time == d, drop=True)
        ds.to_netcdf('./maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/NRT_BFN_'+date.strftime(format)+'.nc', mode = 'w') # Store in daily folder
        ds_renamed = ds.rename({'lon': 'longitude','lat': 'latitude', 'u': 'ugos', 'v': 'vgos'})
        ds_renamed.to_netcdf('./maps_'+name_experiment+'/full_timeseries/BFN_lamta_'+date.strftime(format)+'.nc', mode = 'w') # Store in global product for diagnostics


import sys
from datetime import date
import matplotlib.pyplot as plt

def apply_lamta(name_experiment, currdir, dir_lamta, today, bbox, numdays = 30, delta0 = 0.02, step = 4, bathylvl = -500):
    
    """
    Applies the Lagrangian Analysis of Mixing and Tracer Advection (LAMTA) algorithm to the given experiment.

    Args:
        name_experiment (str): The name of the experiment.
        currdir (str): The current directory.
        dir_lamta (str): The directory of the LAMTA algorithm.
        today (datetime.date): The current date.
        bbox (list): The bounding box coordinates [lon_min, lon_max, lat_min, lat_max].
        numdays (int): The number of days to consider. Default is 30.
        delta0 (float): The initial distance between particles (°). Default is 0.02.
        step (int): The number of integration steps per day. Default is 4.
        bathylvl (int): The bathymetry level from which to study particles resurfacing. Default is -500.

    Returns:
        xr.Dataset: The Lagrangian diagnostic results.

    The function applies the LAMTA algorithm to the given experiment. It loads the required modules and data,
    performs the LAMTA analysis, and saves the Lagrangian diagnostic results. Additionally, it generates and saves
    various plots related to the analysis.
    """

    sys.path.append(dir_lamta)

    from Diagnostics import ParticleSet, Lagrangian
    from Load_nc import load_from_path, loadbathy
    from tools.plot_tools import show_ssh_ug_xinorm, plot_diags_no_sst

    dayv = today.strftime('%Y-%m-%d')
    loni = [bbox[0],bbox[1]]
    lati = [bbox[2],bbox[3]]

    #load bathymetry
    os.chdir(dir_lamta)
    bfield = loadbathy('ETOPO_2022_v1_30s_N90W180_bed.nc',rlon=loni,rlat=lati)
    os.chdir(currdir)

    # Create .nc files with the right names for Lamta
    bfn_output = xr.open_mfdataset(currdir+'/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/NRT_BFN_2*.nc', concat_dim='time', combine='nested')

    # Load field from renamed BFN data in "field"
    first_date = (date.fromisoformat(dayv) - timedelta(days=numdays+3))
    all_days = [(first_date + timedelta(days=x)).strftime('%Y%m%d') for x in range(numdays+4)] # (!) need more dates than numday. e.g. 10 dates for numday = 7
    data_paths = [currdir+'/maps_'+name_experiment+'/full_timeseries/BFN_lamta_'+(first_date + timedelta(days=x)).strftime('%Y%m%d')+'.nc' for x in range(numdays+4)]
    field = load_from_path(data_paths,all_days,unit='deg/d')


    # Trajectories RK4 with u,v
    pset = ParticleSet.from_grid(numdays,loni,lati,delta0,dayv,fieldset=field,mode='backward')
    out = pset.diag(diag=['LLADV','FTLE','OWTRAJ','TIMEFROMBATHY'],method='rk4flat',
                    f=Lagrangian.interpf,numstep=step,
                    coordinates='spherical',numdays=numdays,dayv=dayv,
                    ds=1/6,bathyfield=bfield,
                    bathylvl=bathylvl)

    trjf = out[0] # In case you want to plot particles trajectories - they make pretty rainbow eddies

    diags_results = xr.Dataset(
        data_vars=dict(
            lonf_map=(["lat", "lon"], np.squeeze(out[1]['lonf_map'])),
            latf_map=(["lat", "lon"], np.squeeze(out[1]['latf_map'])),
            ftle=(["lat", "lon"], np.squeeze(out[2]['ftle'])),
            owdisp=(["lat", "lon"], np.squeeze(out[3]['owdisp'])),
            timfb=(["lat", "lon"], np.squeeze(out[4]['timfb'])),
            latfb=(["lat", "lon"], np.squeeze(out[4]['latfb'])),
            lonfb=(["lat", "lon"], np.squeeze(out[4]['lonfb'])),
        ),
        coords=dict(
            lon=(out[0]['lons']),
            lat=(out[0]['lats']),
        ),
        attrs=dict(description="Lagragian diagnostics on provided data.")
    )

    # Save results of lagragian diagnostics
    os.makedirs(currdir+'/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/Lamta/', exist_ok = True)
    diags_results.to_netcdf(currdir+'/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/Lamta/NRT_BFN_lamta_diags'+today.strftime('%Y%m%d')+'.nc', mode = 'w')

    # Make and save Lagrangian diagnostics plots
    fig = show_ssh_ug_xinorm(bfn_output,bbox,np.datetime64(dayv),'BFN output maps')
    fig.savefig(currdir+'/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/Lamta/'+today.strftime('%Y%m%d')+'_BFN_out.png',bbox_inches='tight')
    plt.show

    plot_diags_no_sst(diags_results,bbox,today,bathylvl,save_folder=currdir+'/maps/'+today.strftime('%Y%m%d')+'/Lamta/')
    return diags_results