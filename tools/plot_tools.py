import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
import matplotlib.patches as mpatches
from datetime import timedelta
from pandas import to_datetime
import xarray as xr
import os

def where_is_this(bbox, pad = 4):

    """
    Displays a map showing the specified bounding box.

    Args:
        bbox (list): A list containing the bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat].
        pad (float): Optional. The padding value in degrees to extend the bounding box. Default is 4.

    Returns:
        None

    The function creates a map using Matplotlib and Cartopy to display the specified bounding box. The map is centered
    on the bounding box with a given padding value.
    """

    ax = plt.subplot(projection=ccrs.PlateCarree())
    ax.set_global()
    ax.set_extent([bbox[0]-pad, bbox[1]+pad, bbox[2]-pad, bbox[3]+pad],crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND,color='green',zorder=1)
    ax.add_feature(cfeature.OCEAN,color='lightblue',zorder=0)
    ax.coastlines(lw=0.5)

    ax.add_patch(mpatches.Rectangle(xy=[bbox[0], bbox[2]],
                                    width=bbox[1]-bbox[0],
                                    height=bbox[3]-bbox[2],
                                    edgecolor='red',
                                    facecolor='red',
                                    alpha=0.7,
                                    transform=ccrs.PlateCarree())
                                    )

    plt.show()


# map plotting function to create individual pcolormesh subplot with ax objects
def format_axis_map(data, axes, longitudes, latitudes, bbox = None, colormap='RdBu_r', subplot_title='*insert title here*', cmap_range=None):
    
    """
    Formats and displays a map with data on the specified axes.

    Args:
        data (array-like): The data to be plotted on the map.
        axes (matplotlib.axes.Axes): The axes on which to plot the map.
        longitudes (array-like): The array of longitude values.
        latitudes (array-like): The array of latitude values.
        bbox (list): Optional. The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat]. If not provided, the extent is automatically determined from the longitudes and latitudes. Default is None.
        colormap (str): Optional. The colormap to be used for the data. Default is 'RdBu_r'.
        subplot_title (str): The title of the subplot.
        cmap_range (float or list or tuple): Optional. The range of values for the colormap. If a single float is provided, it is used as both vmin and vmax. If a list or tuple is provided, it should contain the vmin and vmax values. If not provided, the maximum absolute value of the data is used.

    Returns:
        tuple: A tuple containing the formatted axes and the map object.

    The function takes the data, axes, longitudes, and latitudes as inputs and plots the data on the specified
    axes as a map. The map is formatted with the specified colormap, subplot title, and color range. The bounding
    box can be explicitly provided or automatically determined from the longitudes and latitudes. The function
    returns a tuple containing the formatted axes and the map object.
    """
    
    if bbox == None:
        bbox = [np.min(longitudes), np.max(longitudes), np.min(latitudes), np.max(latitudes)]
        
    axes.set_extent([bbox[0], bbox[1], bbox[2], bbox[3]],crs=ccrs.PlateCarree())

    if(cmap_range==None):
        cmap_range = np.nanmax(np.absolute(data))

    if isinstance(cmap_range, (list, tuple, np.ndarray)):
        map = axes.pcolormesh(
                longitudes, 
                latitudes, 
                data, 
                cmap = colormap,
                vmin = cmap_range[0],
                vmax = cmap_range[1]
                )
    else:    
        if(np.nanmin(data) < 0 and np.nanmax(data) > 0):
            map = axes.pcolormesh(
                longitudes, 
                latitudes, 
                data, 
                cmap = colormap,
                vmin = -cmap_range,
                vmax = cmap_range
                )
        else:
            map = axes.pcolormesh(
                longitudes, 
                latitudes, 
                data, 
                cmap = colormap,
                vmax = cmap_range
                )

    axes.add_feature(cfeature.LAND,color='beige',zorder=1)
    axes.add_feature(cfeature.OCEAN,color='lightblue',zorder=0)
    axes.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.1)

    axes.set_xticks(np.arange(bbox[0],bbox[1]+1,2), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axes.xaxis.set_major_formatter(lon_formatter)

    axes.set_yticks(np.arange(bbox[2],bbox[3]+1,2), crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axes.yaxis.set_major_formatter(lat_formatter)

    axes.set_title(subplot_title, loc = 'center')

    return (axes, map)



def format_axis_quiver(data_x, data_y, axes, longitudes, latitudes, bbox = None, subplot_title='*insert title here*'):
    
    """
    Formats and displays a quiver plot on the specified axes.

    Args:
        data_x (array-like): The x-component of the vector data.
        data_y (array-like): The y-component of the vector data.
        axes (matplotlib.axes.Axes): The axes on which to plot the quiver plot.
        longitudes (array-like): The array of longitude values.
        latitudes (array-like): The array of latitude values.
        bbox (list): Optional. The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat]. If not provided, the extent is automatically determined from the longitudes and latitudes. Default is None.
        subplot_title (str): The title of the subplot.

    Returns:
        tuple: A tuple containing the formatted axes and the quiver plot object.

    The function takes the x-component, y-component, axes, longitudes, and latitudes as inputs and plots a quiver
    plot on the specified axes. The quiver plot represents the vector data using arrows. The plot is formatted with
    the specified subplot title and bounding box. The bounding box can be explicitly provided or automatically
    determined from the longitudes and latitudes. The function returns a tuple containing the formatted axes and the
    quiver plot object.
    """
    
    if bbox == None:
        bbox = [np.min(longitudes), np.max(longitudes), np.min(latitudes), np.max(latitudes)]
        
    axes.set_extent([bbox[0], bbox[1], bbox[2], bbox[3]],crs=ccrs.PlateCarree())

    skip = slice(None, None, 2)

    map = axes.quiver(longitudes[(skip)], latitudes[(skip)], data_x[(skip, skip)], data_y[(skip, skip)], scale = 6)

    axes.add_feature(cfeature.LAND,color='beige',zorder=1)
    axes.add_feature(cfeature.OCEAN,color='lightblue',zorder=0)
    axes.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.1)

    axes.set_xticks(np.arange(bbox[0],bbox[1]+1,2), crs=ccrs.PlateCarree())
    lon_formatter = cticker.LongitudeFormatter()
    axes.xaxis.set_major_formatter(lon_formatter)

    axes.set_yticks(np.arange(bbox[2],bbox[3]+1,2), crs=ccrs.PlateCarree())
    lat_formatter = cticker.LatitudeFormatter()
    axes.yaxis.set_major_formatter(lat_formatter)

    axes.set_title(subplot_title, loc = 'center')

    return (axes, map)

def show_ssh_ug_xinorm(ds, bbox, snapshot_time, title = '*insert title here*', s_factor=0.12, title_adjust=1.6):

    """
    Displays a figure with three subplots showing sea surface height, velocity norm, and relative vorticity data.

    Args:
        ds (xarray.Dataset): The dataset containing the required variables.
        bbox (list): The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat].
        snapshot_time (datetime.datetime): The time snapshot to display.
        title (str): The title of the figure.
        s_factor (float): Optional. The scale factor for colorbar size. Default is 0.12.
        title_adjust (float): Optional. The top adjustment value for the figure suptitle. Default is 1.6.

    Returns:
        matplotlib.figure.Figure: The created figure.

    The function takes a dataset, bounding box coordinates, snapshot time, title, and scaling parameters as inputs
    and displays a figure with three subplots. (a) sea surface height (SSH) or absolute dynamic
    topography (ADT) data, (b) velocity norm (Ug) data, (c) realtive vorticity (xi/f).
    Each subplot is formatted with the specified colormap, title, and color range. The function returns the created figure.
    """

    plt.rcParams.update({'font.size': 6})

    if 'ssh' in list(ds.keys()):
        height_var = ds.ssh.sel(time=snapshot_time).data
        height_name = 'SSH'
    elif 'adt' in list(ds.keys()):
        height_var = ds.adt.sel(time=snapshot_time).data
        height_name = 'ADT'
    else:
        return 'Error: no proper surface topography variable'

    u = ds.u.sel(time=snapshot_time).data
    v = ds.v.sel(time=snapshot_time).data
    norm = np.sqrt(u**2 + v**2)
    xi_norm = ds.xi_norm.sel(time=snapshot_time).data
    lons = ds.lon.data
    lats = ds.lat.data

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(height_var, ax1, lons, lats, bbox,'RdBu_r', height_name+' (m)')
    ax1.tick_params('y')
    plt.colorbar(map1, shrink = s_factor, orientation = 'vertical', label = '', ax=ax1)

    ax2, map2 = format_axis_map(norm, ax2, lons, lats, bbox,'viridis', 'Ug (m/s)')
    ax2.tick_params('y', labelleft=False)
    plt.colorbar(map2, shrink = s_factor, orientation = 'vertical', label = '', ax=ax2)

    ax3, map3 = format_axis_map(xi_norm, ax3, lons, lats, bbox,'RdBu_r', r'$\xi$/f', 1) # outlier values at boundaries
    ax3.tick_params('y', labelleft=False)
    plt.colorbar(map3, shrink = s_factor, orientation = 'vertical', label = '', ax=ax3)

    fig.suptitle(title+' on '+np.datetime_as_string(snapshot_time))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust)
    return fig

def show_ug_trio(ds1, ds2, ds3, bbox, snapshot_time, title = '*insert title here*', s_factor=0.12, title_adjust=1.6):

    """
    Displays a figure with three subplots showing velocity norm data from three datasets.

    Args:
        ds1 (xarray.Dataset): The first dataset.
        ds2 (xarray.Dataset): The second dataset.
        ds3 (xarray.Dataset): The third dataset.
        bbox (list): The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat].
        snapshot_time (datetime.datetime): The time snapshot to display.
        title (str): The title of the figure.
        s_factor (float): Optional. The scale factor for colorbar size. Default is 0.12.
        title_adjust (float): Optional. The top adjustment value for the figure suptitle. Default is 1.6.

    Returns:
        None

    The function takes three datasets, a bounding box, snapshot time, title, and scaling parameters as inputs and
    displays a figure with three subplots. Each subplot shows the velocity norm (Ug) data from the respective dataset.
    Each subplot is formatted with the specified colormap, title, and color range.
    """

    plt.rcParams.update({'font.size': 6})

    def get_u_norm(ds):
        norm = np.sqrt(ds.u.sel(time=snapshot_time).data**2 + ds.v.sel(time=snapshot_time).data**2)
        lons = ds.lon.data
        lats = ds.lat.data
        return norm, lons, lats
    
    norm1, lons1, lats1 = get_u_norm(ds1)
    norm2, lons2, lats2 = get_u_norm(ds2)
    norm3, lons3, lats3 = get_u_norm(ds3)

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    def plot_u_norm(norm,ax,lons,lats):
        ax, map = format_axis_map(norm, ax, lons, lats, bbox,'viridis', '')
        ax.tick_params('y', labelleft=False)
        plt.colorbar(map, shrink = s_factor, orientation = 'vertical', label = 'Ug (m/s)', ax=ax)

    plot_u_norm(norm1,ax1,lons1,lats1)
    plot_u_norm(norm2,ax2,lons2,lats2)
    plot_u_norm(norm3,ax3,lons3,lats3)

    fig.suptitle(title+' on '+np.datetime_as_string(snapshot_time))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust)
    plt.show()

def plot_diags_no_sst(diags_results, bbox, today, bathylvl, crop='', s_factor = [0.12,0.12,0.35], title_adjust = [1.6,1.6,1.5], save_folder = './maps/'):
    
    """
    Plots Lamta diagnostic results. 

    Args:
        diags_results: The diagnostic results object (made using Louise Rousslet's algorithm).
        bbox (list): The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat].
        today (datetime.datetime): The current date.
        bathylvl (float): The bathymetry level.
        crop (str): Optional. Crop parameter for the figure filenames. Default is an empty string.
        s_factor (list): Optional. The scale factors for colorbar size. Default is [0.12, 0.12, 0.35].
        title_adjust (list): Optional. The top adjustment values for the figure suptitles. Default is [1.6, 1.6, 1.5].
        save_folder (str): Optional. The folder path to save the figures. Default is './maps/'.

    Returns:
        None

    The function plots Lamta diagnostic results. It generates three sets of subplots for different diagnostic results. 
    Each subplot is formatted with the specified colormap, title, color range, and colorbar label. 
    The figures are saved and displayed but not returned.
    """

    # Creating the FTLE colormap
    cmap = plt.cm.ocean_r
    n=90
    white = np.ones((n,4))
    upper = cmap(np.linspace(0.1, 1, 256-n))
    colors = np.vstack((white, upper))
    ftle_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('fsle_map_white', colors)

    plt.rcParams.update({'font.size': 6})
    title = 'Longitude and latitude advection and Finite-Time Lyapunov Exponents'

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.lonf_map, ax1, diags_results.lon, diags_results.lat, bbox,'BrBG_r', 'Adv lon')
    ax1.tick_params('y')
    ax1.set_aspect('equal')
    plt.colorbar(map1, shrink = s_factor[0], orientation = 'vertical', label = '[°]', ax=ax1)

    ax2, map2 = format_axis_map(diags_results.latf_map, ax2, diags_results.lon, diags_results.lat, bbox,'RdBu_r', 'Adv lat')
    ax2.tick_params('y', labelleft=False)
    ax2.set_aspect('equal')
    plt.colorbar(map2, shrink = s_factor[0], orientation = 'vertical', label = '[°]', ax=ax2)

    ax3, map3 = format_axis_map(diags_results.ftle, ax3, diags_results.lon, diags_results.lat, bbox,ftle_cmap, 'FTLE', [0,0.5])
    ax3.tick_params('y', labelleft=False)
    ax3.set_aspect('equal')
    plt.colorbar(map3, shrink = s_factor[0], orientation = 'vertical', label = '[d-1]', ax=ax3)

    fig.suptitle(title+' on '+today.strftime('%Y%m%d'))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust[0])
    plt.savefig(save_folder+today.strftime('%Y%m%d')+'_lonlatadv_FSLE'+crop+'.png',bbox_inches='tight')
    plt.show()


    title = 'Time, lon, lat from bathy level = '+str(abs(bathylvl))+'m'

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.timfb, ax1, diags_results.lon, diags_results.lat, bbox,'turbo', 'Time from '+str(abs(bathylvl))+'m')
    ax1.tick_params('y')
    ax1.set_aspect('equal')
    plt.colorbar(map1, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax1)

    ax2, map2 = format_axis_map(diags_results.lonfb, ax2, diags_results.lon, diags_results.lat, bbox,'turbo', 'Longitude from '+str(abs(bathylvl))+'m', bbox[0:2])
    ax2.tick_params('y', labelleft=False)
    ax2.set_aspect('equal')
    plt.colorbar(map2, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax2)

    ax3, map3 = format_axis_map(diags_results.latfb, ax3, diags_results.lon, diags_results.lat, bbox,'turbo', 'Latitude from '+str(abs(bathylvl))+'m', bbox[2:])
    ax3.tick_params('y', labelleft=False)
    ax3.set_aspect('equal')
    plt.colorbar(map3, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax3)

    fig.suptitle(title+' on '+today.strftime('%Y%m%d'))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust[1])
    plt.savefig(save_folder+today.strftime('%Y%m%d')+'_diags_bathy'+crop+'.png',bbox_inches='tight')
    plt.show()

    title = 'OW disp'

    fig, ax1 = plt.subplots(nrows=1, ncols=1, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.owdisp, ax1, diags_results.lon, diags_results.lat, bbox,'PuOr', title+' on '+today.strftime('%Y%m%d'), [-15,0])
    plt.colorbar(map1, shrink = s_factor[2], orientation = 'vertical', label = '', ax=ax1)
    ax1.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(save_folder+today.strftime('%Y%m%d')+'_OWdisp'+crop+'.png',bbox_inches='tight')
    plt.show()

def plot_l3_data(bbox, datasets, today, numdays, name_exp):

    """
    Plots the spatial and temporal distribution of L3 data used in the experiment.

    Args:
        bbox (list): The bounding box coordinates in the format [min_lon, max_lon, min_lat, max_lat].
        datasets (list): A list of dataset names.
        today (datetime.datetime): The current date.
        numdays (int): The number of days to consider.
        name_exp (str): The name of the experiment.

    Returns:
        None

    The function plots the spatial and temporal distribution of L3 data for each dataset in the given list. It generates
    two sets of subplots: one for spatial distribution and one for temporal distribution. Each subplot is formatted with
    the specified colormap, color range, and colorbar. The figures are saved and displayed but not returned.
    """

    input_path = './scratch/'+name_exp+'/'

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    #Show tracks in space

    # set number of columns (use 3 to demonstrate the change)
    ncols = 4

    # calculate number of rows
    nrows = len(datasets) // ncols + (len(datasets) % ncols > 0)
    
    plt.figure(figsize=(7*ncols, 7*nrows))
    plt.subplots_adjust(hspace=0.4, wspace=0.7)
    plt.suptitle("Spatial distribution of observations", fontsize=18, y=0.95)

    # loop through the length of tickers and keep track of index
    for n, dataset in enumerate(datasets):
        # add a new subplot iteratively using nrows and cols
        ax = plt.subplot(nrows, ncols, n + 1)

        # filter df and plot ticker on the new subplot axis
        ds = xr.open_mfdataset(input_path+dataset+'*.nc', combine='nested', concat_dim = 'time')
        ds = ds.where((ds['longitude']>bbox[0]) & (ds['longitude']<bbox[1]) & (ds['latitude']>bbox[2]) & (ds['latitude']<bbox[3]),drop = True)
        ds = ds.where((ds['time']>=to_datetime(today-timedelta(days=numdays))) & (ds['time']<=to_datetime(today)), drop = True)
        cmap_range = np.nanmax(np.absolute(ds.SSH))
        a = ax.scatter(ds.SSH.longitude, ds.SSH.latitude, c=ds.SSH, cmap = 'RdBu_r', vmin = -cmap_range, vmax = cmap_range)

        # chart formatting
        ax.set_title(dataset)
        ax.set_aspect('equal', 'box')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(a, cax=cax)

    plt.tight_layout()
    os.makedirs('./maps_'+name_exp+'/'+today.strftime('%Y%m%d')+'/', exist_ok = True)
    plt.savefig('./maps_'+name_exp+'/'+today.strftime('%Y%m%d')+'/L3_data_map'+'.png',bbox_inches='tight')
    plt.show()

    #Show tracks in time
    plt.figure(figsize=(7*ncols, 7*nrows))
    plt.subplots_adjust(hspace=0.4, wspace=0.5)
    plt.suptitle("Temporal distribution of observations", fontsize=18, y=0.95)

    # loop through the length of tickers and keep track of index
    for n, dataset in enumerate(datasets):
        # add a new subplot iteratively using nrows and cols
        ax = plt.subplot(nrows, ncols, n + 1)

        # filter df and plot ticker on the new subplot axis
        ds = xr.open_mfdataset(input_path+dataset+'*.nc', combine='nested', concat_dim = 'time')
        ds = ds.where((ds['longitude']>bbox[0]) & (ds['longitude']<bbox[1]) & (ds['latitude']>bbox[2]) & (ds['latitude']<bbox[3]),drop = True)
        ds = ds.where((ds['time']>=to_datetime(today-timedelta(days=numdays))) & (ds['time']<=to_datetime(today)), drop = True)
        ax.scatter(ds.time,ds.SSH)

        # chart formatting
        ax.set_title(dataset)
        plt.xticks(rotation=45)

    plt.tight_layout()
    plt.savefig('./maps_'+name_exp+'/'+today.strftime('%Y%m%d')+'/L3_data_timeseries'+'.png',bbox_inches='tight')
    plt.show()

import glob

def plot_alongtrack_rmse(input_path, name_exp, date_folder):

    """
    Plots the RMSE between maps and observation tracks.

    Args:
        input_path (str): The path to the observation tracks.
        name_exp (str): The name of the experiment.
        date_folder (str): The date folder.

    Returns:
        None

    The function opens the necessary map datasets and observation tracks. It calculates the RMSE between the maps and
    observation tracks and plots the results. The figures are displayed but not returned.
    """

    bfn_map = xr.open_mfdataset('./output_'+name_exp+'/'+date_folder+'/*.nc', combine='nested', concat_dim = 'time')
    duacs_map = xr.open_mfdataset('./input_'+name_exp+'/'+date_folder+'/dataset-duacs-nrt-global-merged-allsat-phy-l4/*.nc', combine='nested', concat_dim = 'time')

    tracks = glob.glob(input_path+'obs*')
    bfn_rmse = np.zeros(len(tracks))
    duacs_rmse = np.zeros(len(tracks))

    def rmse(predictions, targets):
        return np.sqrt(np.mean((predictions-targets)**2))

    for n, track in enumerate(tracks) :
        obs = xr.open_dataset(track)
        bfn_on_track = bfn_map['ssh'].interp(lon = obs.longitude, lat = obs.latitude, time = obs.time)
        duacs_on_track = duacs_map['adt'].interp(longitude = obs.longitude, latitude = obs.latitude, time = obs.time)

        bfn_rmse[n]=rmse(bfn_on_track,obs['SSH'])
        duacs_rmse[n]=rmse(duacs_on_track,obs['SSH'])

    plt.plot(duacs_rmse, label='duacs (total : '+str(np.nansum(duacs_rmse))+')')
    plt.plot(bfn_rmse, label='bfn (total : '+str(np.nansum(bfn_rmse))+')')
    plt.title('RMSE between maps and observation tracks')
    plt.legend()
    plt.savefig('./maps_'+name_exp+'/'+date_folder+'/RMSE_to_obs_tracks.png',bbox_inches='tight')
    plt.show()

def plot_25_random_tracks(input_path, name_exp, date_folder):

    """
    Plots 25 randomly chosen observation tracks along with map data.

    Args:
        input_path (str): The path to the observation tracks.
        name_exp (str): The name of the experiment.
        date_folder (str): The date folder.

    Returns:
        None

    The function opens the necessary map datasets and selects 25 random observation tracks. It plots each observation
    track along with the corresponding map data. The figures are displayed but not returned.
    """

    plt.figure(figsize=(15, 15))
    plt.subplots_adjust(hspace=0.4)

    tracks = np.random.choice(glob.glob(input_path+'obs*'), 25)
    # To choose tracks explicitely, use :
    #tracks = glob.glob(input_path+'obs*')[0:25]

    bfn_map = xr.open_mfdataset('./output_'+name_exp+'/'+date_folder+'/*.nc', combine='nested', concat_dim = 'time')
    duacs_map = xr.open_mfdataset('./input_'+name_exp+'/'+date_folder+'/dataset-duacs-nrt-global-merged-allsat-phy-l4/*.nc', combine='nested', concat_dim = 'time')

    # set number of columns and rows
    ncols = 5
    nrows = len(tracks) // ncols + (len(tracks) % ncols > 0)

    for n, track in enumerate(tracks) :
        # add a new subplot iteratively using nrows and cols
        ax = plt.subplot(nrows, ncols, n + 1)

        obs = xr.open_dataset(track)
        bfn_on_track = bfn_map['ssh'].interp(lon = obs.longitude, lat = obs.latitude, time = obs.time)
        duacs_on_track = duacs_map['adt'].interp(longitude = obs.longitude, latitude = obs.latitude, time = obs.time)
        
        ax.plot(obs['SSH'].time.values, obs['SSH'].values, label = 'alongtrack observation')
        ax.plot(bfn_on_track.time.values, bfn_on_track.values, label = 'BFN ouput interpolated alongtrack')
        ax.plot(duacs_on_track.time.values, duacs_on_track.values, label = 'DUACS L4 interpolated alongtrack')
        plt.xticks(rotation=45)

    plt.legend()
    plt.savefig('./maps_'+name_exp+'/'+date_folder+'/25_random_tracks_check.png',bbox_inches='tight')
    plt.show()

def make_comp_plot(name_exp, date_folder, bbox, snap_time, bfn_map = None, duacs_map_crop = None):

    """
    Creates a comparison plot of DUACS and BFN maps for a given snapshot time.

    Args:
        name_exp (str): The name of the experiment.
        date_folder (str): The date folder.
        bbox (list): The bounding box coordinates [lon_min, lon_max, lat_min, lat_max].
        snap_time (datetime): The snapshot time.
        bfn_map (xr.Dataset, optional): The BFN map dataset. If not provided, it is loaded from a default path.
        duacs_map_crop (xr.Dataset, optional): The cropped DUACS map dataset. If not provided, it is loaded from a default path.

    Returns:
        None

    The function creates a comparison plot of DUACS and BFN maps for the specified snapshot time. It opens the necessary
    datasets, sets up the figure and subplots, and formats the axis maps. The resulting plot is displayed.
    """

    if bfn_map is None:
        bfn_map = xr.open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT-BFN-v2/maps_'+name_exp+'/'+date_folder+'/*.nc', combine='nested', concat_dim = 'time')
        duacs_map = xr.open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT-BFN-v2/input_'+name_exp+'/'+date_folder+'/dataset-duacs-nrt-global-merged-allsat-phy-l4/*.nc', combine='nested', concat_dim = 'time')
        duacs_map_crop = duacs_map.where((duacs_map['longitude']>bbox[0]) & (duacs_map['longitude']<bbox[1]) & (duacs_map['latitude']>bbox[2]) & (duacs_map['latitude']<bbox[3]),drop = True)
    
    duacs = duacs_map_crop.sel(time=to_datetime(snap_time), method = 'nearest')
    bfn = bfn_map.sel(time=to_datetime(snap_time), method = 'nearest')
    
    fig, ((ax1, ax2),(ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    cmap_range_ssh = np.nanmax([np.nanmax(abs(duacs.adt.values)), np.nanmax(abs(bfn.ssh.values))])

    ax1, map1 = format_axis_map(duacs.adt, ax1, duacs.longitude, duacs.latitude, colormap='RdBu_r', subplot_title='DUACS', cmap_range=[-cmap_range_ssh,cmap_range_ssh])
    ax1.tick_params('x', labelbottom=False)
    plt.colorbar(map1, shrink = 0.4, orientation = 'vertical', label = 'ADT (m)', ax=ax1)

    ax2, map2 = format_axis_map(bfn.ssh, ax2, bfn.lon, bfn.lat, colormap='RdBu_r', subplot_title='BFN-QG', cmap_range=[-cmap_range_ssh,cmap_range_ssh])
    ax2.tick_params('y', labelleft=False)
    ax2.tick_params('x', labelbottom=False)
    plt.colorbar(map2, shrink = 0.4, orientation = 'vertical', label = 'ADT (m)', ax=ax2)

    cmap_range_uv = np.nanmax([np.nanmax(np.sqrt(duacs.ugos.values**2+duacs.vgos.values**2)), np.nanmax(np.sqrt(bfn.u.values**2+bfn.v.values**2))])

    ax3, map3 = format_axis_map(np.sqrt(duacs.ugos**2+duacs.vgos**2), ax3, duacs.longitude, duacs.latitude, colormap='viridis', subplot_title='', cmap_range=[0,cmap_range_uv])
    ax3.tick_params('x', labelrotation=45)
    plt.colorbar(map3, shrink = 0.4, orientation = 'vertical', label = 'sqrt(u²+v²) (m/s)', ax=ax3)
    

    ax4, map4 = format_axis_map(np.sqrt(bfn.u**2+bfn.v**2), ax4, bfn.lon, bfn.lat, colormap='viridis', subplot_title='', cmap_range=[0,cmap_range_uv])
    ax4.tick_params('y', labelleft=False)
    ax4.tick_params('x', labelrotation=45)
    plt.colorbar(map4, shrink = 0.4, orientation = 'vertical', label = 'sqrt(u²+v²) (m/s)', ax=ax4)

    fig.suptitle('Surface elevation & norm of geostrophic velocity, '+snap_time.strftime('%Y-%m-%d'))
    fig.tight_layout()

def plot_duacs_comp (init_date, name_experiment, today, bbox, make_duacs_comp = 'today'):

    """
    Calls the function plotting the comparison between DUACS and BFN maps depending on the specified option.

    Args:
        init_date (datetime.date): The initial date.
        name_experiment (str): The name of the experiment.
        today (datetime.date): The current date.
        bbox (list): The bounding box coordinates [lon_min, lon_max, lat_min, lat_max].
        make_duacs_comp (str, optional): The type of DUACS comparison to make. Default is 'today'.

    Returns:
        None

    The function plots the comparison between DUACS and BFN maps based on the specified options. If 'make_duacs_comp'
    is set to 'today', it calls the 'make_comp_plot' function for the current date. If set to 'interactive', it allows
    for interactive selection of a specific date to compare. If set to 'none', no DUACS comparison is made. Otherwise,
    if a specific date is provided in 'YYYY-MM-DD' format, it calls the 'make_comp_plot' function for that date.
    """

    if make_duacs_comp == 'today':
        from tools.plot_tools import make_comp_plot
        make_comp_plot(name_experiment, today.strftime('%Y%m%d'), bbox, today)
        plt.savefig('./maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/DUACS_L4_comp_'+today.strftime('%Y%m%d')+'.png',bbox_inches='tight')

    elif make_duacs_comp == 'interactive':
        from tools.plot_tools import make_comp_plot
        import ipywidgets
        from numpy import arange
        from pandas import to_datetime
        from xarray import open_mfdataset

        dates_series = to_datetime(arange(init_date,today+timedelta(days=1),timedelta(days=1)))
        bfn_map = open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT-BFN-v2/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/*.nc', combine='nested', concat_dim = 'time')
        duacs_map = open_mfdataset('/bettik/PROJECTS/pr-data-ocean/stellaa/NRT-BFN-v2/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/dataset-duacs-nrt-global-merged-allsat-phy-l4/*.nc', combine='nested', concat_dim = 'time')
        duacs_map_crop = duacs_map.where((duacs_map['longitude']>bbox[0]) & (duacs_map['longitude']<bbox[1]) & (duacs_map['latitude']>bbox[2]) & (duacs_map['latitude']<bbox[3]),drop = True)

        def plot_interactive_comp(i):
            make_comp_plot(name_exp=name_experiment, date_folder=today.strftime('%Y%m%d'), bbox=bbox, snap_time=dates_series[i], bfn_map=bfn_map, duacs_map_crop=duacs_map_crop)

        w = ipywidgets.interactive(plot_interactive_comp, i=(0, len(dates_series)-1))
        display(w)
    elif make_duacs_comp == 'none':
        print('No duacs comparaison requested')

    else:
        try:
            from tools.plot_tools import make_comp_plot
            from datetime import datetime, date, time
            plot_date = datetime.combine(date.fromisoformat(make_duacs_comp), time())
            make_comp_plot(name_experiment, today.strftime('%Y%m%d'), bbox, plot_date)
            plt.savefig('./maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/DUACS_L4_comp_'+plot_date.strftime('%Y%m%d')+'.png',bbox_inches='tight')
        except ValueError:
            print("invalid value for make_duacs_comp. Choose from 'today', 'interactive', 'none', or a date in 'YYYY-MM-DD' format.")