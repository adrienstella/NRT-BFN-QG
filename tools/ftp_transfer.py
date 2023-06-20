def ftp_cmems_download_month(ftp, product_path, dataset, year, month):

    """ Download all data files in a month repository on an ftp server

    Args:
        ftp (ftp object): FTP connexion
        product_path (string): Name of product to download
        dataset (string): Name of dataset to download within product
        year (string): year of data to download
        month (string): month of data to download

    Returns:
        None
    """

    # Change CMEMS directory to the product we want (find the path in the catalogue)
    ftp.cwd(product_path+dataset)
    ftp.cwd(str(year))
    
    if(len(month)<2):
        month = '0'+month

    if month in ftp.nlst():    
        ftp.cwd(month)

        # Set the name of the file to download
        filenames = ftp.nlst()
        for filename in filenames:
            print('Retreiving data for '+filename)

            # Download the file
            ftp.retrbinary("RETR "+filename, open(filename, 'wb').write)
    
    else:
        print('[Warning] No data folder for month '+ month + ' yet.')


import os
from ftplib import error_perm

def place_files(ftp, path):
    # A function that uploads everything in path (subfolders included) to the directory to which ftp is connected.
    for name in os.listdir(path):
        localpath = os.path.join(path, name)
        if os.path.isfile(localpath):
            print("STOR", name, localpath)
            ftp.storbinary('STOR ' + name, open(localpath,'rb'))
        elif os.path.isdir(localpath):
            print("MKD", name)

            try:
                ftp.mkd(name)

            # ignore "directory already exists"
            except error_perm as e:
                if not e.args[0].startswith('550'): 
                    raise

            print("CWD", name)
            ftp.cwd(name)
            place_files(ftp, localpath)           
            print("CWD", "..")
            ftp.cwd("..")

import secretcodes
from ftplib import FTP
def ftp_to_ifremer(today, currdir):

    # Set user name and password
    username = secretcodes.ifremer_username
    password = secretcodes.ifremer_password

    # Connect to the ftp server
    ftp = FTP('ftp.ifremer.fr',username,password)
    ftp.cwd('MEDSSH_BFN')

    ftp.mkd(today.strftime('%Y%m%d'))
    ftp.cwd(today.strftime('%Y%m%d'))

    repo_to_upload = currdir+'/maps/'+today.strftime('%Y%m%d')+'/'
    place_files(ftp, repo_to_upload)

    print('Closing connexion :')
    ftp.quit()

import numpy as np
from datetime import timedelta
def download_nadirs_cmems(currdir, today, numdays, datasets, dataset_l4):
    # Set user name and password
    username = secretcodes.cmems_username
    password = secretcodes.cmems_password

    # Connect to the ftp server
    ftp = FTP('nrt.cmems-du.eu',username,password)

    # Choose dates to download
    first_day = today - timedelta(days=numdays)
    # Download all L3 products

    for i in np.arange(0,len(datasets)):

        os.makedirs(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+datasets[i], exist_ok = True)
        os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+datasets[i])

        ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L3_NRT_OBSERVATIONS_008_044/', datasets[i], str(today.year), str(today.month))
        os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+datasets[i])

        if (today.month != first_day.month): # (!) only works when the assimilation duration doesn't span over more than two months
                ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L3_NRT_OBSERVATIONS_008_044/', datasets[i], str(first_day.year), str(first_day.month))
                os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+datasets[i])

    os.chdir(currdir)
    print('Obs data downloaded successfully')

    # Download DUACS L4 product
    print('Retreiving data for dataset '+dataset_l4)

    os.makedirs(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+dataset_l4, exist_ok = True)
    os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+dataset_l4)

    ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/', dataset_l4, str(today.year), str(today.month))
    os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+dataset_l4)

    if (today.month != first_day.month): # (!) only works when the assimilation duration doesn't span over more than two months
        ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/', dataset_l4, str(first_day.year), str(first_day.month))

    os.chdir(currdir)
    print('DUACS L4 data downloaded successfully')

    ftp.quit()

def download_mdt(currdir, dataset_mdt):
    # Set user name and password
    username = secretcodes.cmems_username
    password = secretcodes.cmems_password

    # Connect to the ftp server
    ftp = FTP('my.cmems-du.eu',username,password)

    # Download DUACS L4 product
    print('Downloading MDT')

    os.makedirs(currdir+'/input/', exist_ok = True)
    os.chdir(currdir+'/input/')

    ftp.retrbinary("RETR "+'/Core/SEALEVEL_GLO_PHY_MDT_008_063/cnes_obs-sl_glo_phy-mdt_my_0.125deg_P20Y/'+dataset_mdt, open(dataset_mdt, 'wb').write)

    os.chdir(currdir)
    print('MDT downloaded successfully')

    ftp.quit()

import re
def download_swot_nadir(currdir, today):
    # Download SWOT nadir L3 data from AVISO
    # Set user name and password
    username = secretcodes.swot_username
    password = secretcodes.swot_password

    # Connect to the ftp server
    ftp = FTP('ftp-access.aviso.altimetry.fr',username,password)
    ftp.cwd('/data/Data/ALTI/DUACS_SWOT_Nadir/L3_Along_track')
    filenames = ftp.nlst()

    # Download SWOT nadir product
    dataset_swot_n = 'nrt_global_swonc_phy_l3_1hz'
    print('Retreiving data for dataset '+dataset_swot_n)

    os.makedirs(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+dataset_swot_n, exist_ok = True)
    os.chdir(currdir+'/input/'+today.strftime('%Y%m%d')+'/'+dataset_swot_n)

    # Set the name of the file to download
    for filename in filenames:
        print('Retreiving data for '+filename)

        # Download the file
        ftp.retrbinary("RETR "+filename, open(filename, 'wb').write)

    print(dataset_swot_n+' data downloaded successfully')
    ftp.quit()
    
    # Sort SWOT nadir files to keep only the right ones
    # Get all combinations of data dates and upload dates from filenames
    i=0
    all_date_combos = np.zeros((len(filenames),2))
    for filename in filenames:
        all_date_combos[i,0] = re.findall(r'\d+', filename)[2]
        all_date_combos[i,1] = re.findall(r'\d+', filename)[3]
        i=i+1
    all_date_combos = all_date_combos[all_date_combos[:,0].argsort()]

    # Makes a list of all data dates available, regardless of upload date
    dates_data = np.zeros((0,0))
    for date in all_date_combos[:,0]:
        if not(any(dates_data==date)):
            dates_data = np.append(dates_data, date)

    # For each data date, choose the latest upload date
    i=0
    most_recent_uploads = np.zeros((len(dates_data),2))
    for date in dates_data:
        date_i_uploads = all_date_combos[all_date_combos[:,0]==date, :]
        most_recent_uploads[i,0] = date 
        most_recent_uploads[i,1] = max(date_i_uploads[:,1])
        # print('At date '+str(most_recent_uploads[i,0])+', uploads available : '+str(date_i_uploads[:,1])+', keeping : '+str(most_recent_uploads[i,1]))
        i=i+1
    # most_recent_uploads # this contains pairs of data dates and upload dates to KEEP. all other should be deleted.

    keeper_names = []
    for i in np.arange(0,len(most_recent_uploads)):
        keeper_names = np.append(keeper_names, 'nrt_global_swonc_phy_l3_1hz_'+str(int(most_recent_uploads[i,0]))+'_'+str(int(most_recent_uploads[i,1]))+'.nc')

    for filename in filenames:
        if not(filename in keeper_names):
            print('[deleting] '+filename)
            os.remove(filename)
        else:
            print('[keeping] '+filename)

    print('[SWOT nadir input files ready]')
    os.chdir(currdir)