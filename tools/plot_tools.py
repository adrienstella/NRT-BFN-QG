import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
import matplotlib.patches as mpatches

def where_is_this(bbox, pad = 4):
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



# map plotting function to create individual vector subplot with ax objects
def format_axis_quiver(data_x, data_y, axes, longitudes, latitudes, bbox = None, subplot_title='*insert title here*'):
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


def plot_all_diags(diags_results, bbox, today, crop='', s_factor = [0.12,0.12,0.22], title_adjust = [1.6,1.6,1.5], bathylvl = 700):
    
    plt.rcParams.update({'font.size': 6})
    title = 'Longitude and latitude advection and Finite-Time Lyapunov Exponents'

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.lonf_map, ax1, diags_results.lon, diags_results.lat, bbox,'BrBG_r', 'Adv lon')
    ax1.tick_params('y')
    plt.colorbar(map1, shrink = s_factor[0], orientation = 'vertical', label = '', ax=ax1)

    ax2, map2 = format_axis_map(diags_results.latf_map, ax2, diags_results.lon, diags_results.lat, bbox,'RdBu_r', 'Adv lat')
    ax2.tick_params('y', labelleft=False)
    plt.colorbar(map2, shrink = s_factor[0], orientation = 'vertical', label = '', ax=ax2)

    ax3, map3 = format_axis_map(diags_results.ftle, ax3, diags_results.lon, diags_results.lat, bbox,'inferno', 'FTLE', [0,1])
    ax3.tick_params('y', labelleft=False)
    plt.colorbar(map3, shrink = s_factor[0], orientation = 'vertical', label = '', ax=ax3)

    fig.suptitle(title+' on '+today.strftime('%Y%m%d'))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust[0])
    plt.savefig('./maps/'+today.strftime('%Y%m%d')+'/'+today.strftime('%Y%m%d')+'_lonlatadv_FSLE'+crop+'.png',bbox_inches='tight')
    plt.show()


    title = 'Time, lon, lat from bathy level = '+str(abs(bathylvl))+'m'

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.timfb, ax1, diags_results.lon, diags_results.lat, bbox,'turbo', 'Time from '+str(abs(bathylvl))+'m')
    ax1.tick_params('y')
    plt.colorbar(map1, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax1)

    ax2, map2 = format_axis_map(diags_results.lonfb, ax2, diags_results.lon, diags_results.lat, bbox,'turbo', 'Longitude from '+str(abs(bathylvl))+'m', bbox[0:2])
    ax2.tick_params('y', labelleft=False)
    plt.colorbar(map2, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax2)

    ax3, map3 = format_axis_map(diags_results.latfb, ax3, diags_results.lon, diags_results.lat, bbox,'turbo', 'Latitude from '+str(abs(bathylvl))+'m', bbox[2:])
    ax3.tick_params('y', labelleft=False)
    plt.colorbar(map3, shrink = s_factor[1], orientation = 'vertical', label = '', ax=ax3)

    fig.suptitle(title+' on '+today.strftime('%Y%m%d'))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust[1])
    plt.savefig('./maps/'+today.strftime('%Y%m%d')+'/'+today.strftime('%Y%m%d')+'_diags_bathy'+crop+'.png',bbox_inches='tight')
    plt.show()

    title = 'OW disp & SST advection'

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, subplot_kw={'projection': ccrs.PlateCarree()},dpi=300)

    ax1, map1 = format_axis_map(diags_results.owdisp, ax1, diags_results.lon, diags_results.lat, bbox,'PuOr', 'OWdisp', [-15,0])
    ax1.tick_params('y')
    plt.colorbar(map1, shrink = s_factor[2], orientation = 'vertical', label = '', ax=ax1)

    ax2, map2 = format_axis_map(diags_results.sstadv, ax2, diags_results.lon, diags_results.lat, bbox,'inferno', 'SST advection', [14,20])
    ax2.tick_params('y', labelleft=False)
    plt.colorbar(map2, shrink = s_factor[2], orientation = 'vertical', label = '', ax=ax2)

    fig.suptitle(title+' on '+today.strftime('%Y%m%d'))
    plt.tight_layout()
    fig.subplots_adjust(top=title_adjust[2])
    plt.savefig('./maps/'+today.strftime('%Y%m%d')+'/'+today.strftime('%Y%m%d')+'_OWdisp_SST'+crop+'.png',bbox_inches='tight')
    plt.show()

def plot_diags_no_sst(diags_results, bbox, today, bathylvl, crop='', s_factor = [0.12,0.12,0.35], title_adjust = [1.6,1.6,1.5], save_folder = './maps/'):
    
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