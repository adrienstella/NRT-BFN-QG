import xarray as xr
import numpy as np
import pyinterp.backends.xarray
import pyinterp.fill
from pandas import to_datetime
from datetime import timedelta
import os

def compute_filled_map(BC_data_path, save_to, bbox):
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


def nc_processing(today):

    # A few physical parameters for computations
    g=9.81 # m/s2
    earth_rad=6.371e6 # 6371km
    earth_w=2*np.pi/86400 # rad/s

    # Load the last week and take daily averages
    bfn_output = xr.open_mfdataset('./output/'+today.strftime('%Y%m%d')+'/*.nc', concat_dim='time', combine='nested')
    bfn_output_dates = bfn_output.assign(time=to_datetime(bfn_output.time.dt.date))
    last_week = bfn_output_dates.where(bfn_output_dates.time >= to_datetime(today-timedelta(days=6)), drop =True)
    daily_mean_ssh = last_week.groupby("time").mean("time")

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

    # Save final netcdf files - this can be adapted depending on the needs of the receiver (variable names, how many files, etc.)
    os.makedirs('./maps/'+today.strftime('%Y%m%d')+'/', exist_ok = True)
    os.makedirs('./maps/full_timeseries/', exist_ok = True)
    for d in daily_mean_ssh.time.values:
        date = to_datetime(d)
        ds=daily_mean_ssh.where(daily_mean_ssh.time == d, drop=True)
        ds.to_netcdf('./maps/'+today.strftime('%Y%m%d')+'/NRT_BFN_'+date.strftime('%Y%m%d')+'.nc', mode = 'w') # Store in daily folder
        ds_renamed = ds.rename({'lon': 'longitude','lat': 'latitude', 'u': 'ugos', 'v': 'vgos'})
        ds_renamed.to_netcdf('./maps/full_timeseries/BFN_lamta_'+date.strftime('%Y%m%d')+'.nc', mode = 'w') # Store in global product for diagnostics


import sys
from datetime import date
import matplotlib.pyplot as plt
def apply_lamta(currdir, dir_lamta, today, bbox, numdays = 30, delta0 = 0.02, step = 4, bathylvl = -500):
    
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
    bfn_output = xr.open_mfdataset(currdir+'/maps/'+today.strftime('%Y%m%d')+'/NRT_BFN_2*.nc', concat_dim='time', combine='nested')

    # Load field from renamed BFN data in "field"
    first_date = (date.fromisoformat(dayv) - timedelta(days=numdays+3))
    all_days = [(first_date + timedelta(days=x)).strftime('%Y%m%d') for x in range(numdays+4)] # (!) need more dates than numday. e.g. 10 dates for numday = 7
    data_paths = [currdir+'/maps/full_timeseries/BFN_lamta_'+(first_date + timedelta(days=x)).strftime('%Y%m%d')+'.nc' for x in range(numdays+4)]
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
    os.makedirs(currdir+'/maps/'+today.strftime('%Y%m%d')+'/Lamta/', exist_ok = True)
    diags_results.to_netcdf(currdir+'/maps/'+today.strftime('%Y%m%d')+'/Lamta/NRT_BFN_lamta_diags'+today.strftime('%Y%m%d')+'.nc', mode = 'w')

    # Make and save Lagrangian diagnostics plots
    fig = show_ssh_ug_xinorm(bfn_output,bbox,np.datetime64(dayv),'BFN output maps')
    fig.savefig(currdir+'/maps/'+today.strftime('%Y%m%d')+'/Lamta/'+today.strftime('%Y%m%d')+'_BFN_out.png',bbox_inches='tight')
    plt.show

    plot_diags_no_sst(diags_results,bbox,today,bathylvl,save_folder=currdir+'/maps/'+today.strftime('%Y%m%d')+'/Lamta/')
    return diags_results

# Check MDT is ready
def make_mdt(currdir, bbox, dataset_mdt = 'mdt_hybrid_cnes_cls18_cmems2020_global.nc'):
    if(not(os.path.isfile(currdir+'/input/cnes_mdt_local.nc'))):
        from tools.ftp_transfer import download_mdt
        download_mdt(currdir, dataset_mdt)
        ds = xr.open_dataset(currdir+'/input/'+dataset_mdt)
        ds = ds.sel(longitude = slice(bbox[0],bbox[1]), latitude = slice(bbox[2],bbox[3]))
        mdt = ds.mdt
        mdt.to_netcdf(currdir+'/input/cnes_mdt_local.nc')