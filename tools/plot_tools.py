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
    plt.show()