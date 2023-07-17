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
    if(year in ftp.nlst()):
        ftp.cwd(year)
    
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
    else:
        print('[Warning] No data folder for year '+ year + ' yet.')


import os
from ftplib import error_perm

def place_files(ftp, path):
    """
    Uploads all files and subfolders in the given path to the directory connected to the FTP server.

    Args:
        ftp (ftplib.FTP): An FTP object, already connected to the target directory.
        path (str): The local path to the directory containing files and subfolders to upload.

    Returns:
        None

    Raises:
        ftplib.error_perm: If there is an error during FTP operations.

    The function iterates over the files and subfolders in the specified local path. For each file,
    it uploads the file to the connected FTP directory using the STOR command. For each subfolder,
    it creates a corresponding directory on the FTP server using the MKD command. The function then
    recursively calls itself for each subfolder to upload its contents.
    """
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
def ftp_to_ifremer(name_experiment, today, currdir):

    """
    Uploads files to the IFREMER FTP server for a specific experiment.

    Args:
        name_experiment (str): The name of the experiment.
        today (datetime.datetime): The current experiment date.
        currdir (str): The current working directory.

    Returns:
        None

    The function connects to the IFREMER FTP server using the provided credentials and navigates to
    the 'MEDSSH_BFN' directory. It creates a subdirectory with the current date in the 'MEDSSH_BFN'
    directory, then, it calls the 'place_files' function to upload the files to the FTP server. 
    Finally, it closes the FTP connection.

    Note:
        - The 'secretcodes' module must contain the 'ifremer_username' and 'ifremer_password' variables
          with the appropriate FTP credentials.
    """

    # Set user name and password
    username = secretcodes.ifremer_username
    password = secretcodes.ifremer_password

    # Connect to the ftp server
    ftp = FTP('ftp.ifremer.fr',username,password)
    ftp.cwd('MEDSSH_BFN')

    ftp.mkd(today.strftime('%Y%m%d'))
    ftp.cwd(today.strftime('%Y%m%d'))

    repo_to_upload = currdir+'/maps_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'
    place_files(ftp, repo_to_upload)

    print('Closing connexion :')
    ftp.quit()

import numpy as np
from datetime import timedelta

def download_nadirs_cmems(name_experiment, currdir, today, numdays, datasets, dataset_l4):

    """
    Downloads nadir data from the CMEMS FTP server for a specific experiment.

    Args:
        name_experiment (str): The name of the experiment.
        currdir (str): The current working directory.
        today (datetime.datetime): The current experiment date.
        numdays (int): The minimum number of days to go back for data download.
        datasets (list): A list of dataset names to download (alongtrack data).
        dataset_l4 (str): The name of the L4 dataset to download.

    Returns:
        None

    The function connects to the CMEMS FTP server using the provided credentials and downloads nadir
    data for the specified datasets. It creates the necessary directories to store the downloaded data
    within the experiment's input folder. It then iterates over the months starting from the
    current month and goes back as many months as needed to cover 'numdays'. 

    Note:
        - The 'secretcodes' module must contain the 'cmems_username' and 'cmems_password' variables
          with the appropriate FTP credentials.
    """

    # Set user name and password
    username = secretcodes.cmems_username
    password = secretcodes.cmems_password

    # Connect to the ftp server
    ftp = FTP('nrt.cmems-du.eu',username,password)

    # Choose dates to download
    first_day = today - timedelta(days=numdays)
    # Download all L3 products

    for i in np.arange(0,len(datasets)):

        os.makedirs(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+datasets[i], exist_ok = True)
        os.chdir(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+datasets[i])

        # go through each month folder in case obs span over several months
        m=0
        month = today.month
        year = today.year
        while (((month != first_day.month) | (year != first_day.year) | (m==0))):
            month = (today.month-m)%12
            if month == 0:
                month = 12
                year = year-1
            ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L3_NRT_OBSERVATIONS_008_044/', datasets[i], str(year), str(month))
            os.chdir(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+datasets[i])
            m += 1

    os.chdir(currdir)
    print('Obs data downloaded successfully')

    # Download DUACS L4 product
    print('Retreiving data for dataset '+dataset_l4)

    os.makedirs(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+dataset_l4, exist_ok = True)
    os.chdir(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+dataset_l4)
    
    m=0
    month = today.month
    year = today.year
    while (((month != first_day.month) | (year != first_day.year) | (m==0))):
        month = (today.month-m)%12
        if month == 0:
            month = 12
            year = year-1
        ftp_cmems_download_month(ftp, '/Core/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/', dataset_l4, str(year), str(month))
        os.chdir(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+dataset_l4)
        m += 1

    os.chdir(currdir)
    print('DUACS L4 data downloaded successfully')

    ftp.quit()

def download_mdt(name_experiment, currdir, dataset_mdt):

    """
    Downloads Mean Dynamic Topography (MDT) data from the CMEMS FTP server for a specific experiment.

    Args:
        name_experiment (str): The name of the experiment.
        currdir (str): The current working directory.
        dataset_mdt (str): The name of the MDT dataset to download.

    Returns:
        None

    The function connects to the CMEMS FTP server using the provided credentials and downloads the Mean
    Dynamic Topography (MDT) dataset to the specified experiment's input folder. 

    Note:
        - The 'secretcodes' module must contain the 'cmems_username' and 'cmems_password' variables
          with the appropriate FTP credentials.
    """

    # Set user name and password
    username = secretcodes.cmems_username
    password = secretcodes.cmems_password

    # Connect to the ftp server
    ftp = FTP('my.cmems-du.eu',username,password)

    # Download DUACS L4 product
    print('Downloading MDT')

    os.makedirs(currdir+'/input_'+name_experiment+'/', exist_ok = True)
    os.chdir(currdir+'/input_'+name_experiment+'/')

    ftp.retrbinary("RETR "+"/Core/SEALEVEL_GLO_PHY_MDT_008_063/cnes_obs-sl_glo_phy-mdt_my_0.125deg_P20Y/"+dataset_mdt, open(dataset_mdt, 'wb').write)

    os.chdir(currdir)
    print('MDT downloaded successfully')

    ftp.quit()

import re

def download_swot_nadir(name_experiment, currdir, today):

    """
    Downloads SWOT nadir L3 data from AVISO for a specific experiment.

    Args:
        name_experiment (str): The name of the experiment.
        currdir (str): The current working directory.
        today (datetime.datetime): The current date.

    Returns:
        None

    The function connects to the AVISO FTP server using the provided credentials and downloads SWOT nadir
    L3 data for the specified experiment. It retrieves the list of available filenames from the FTP server.
    It iterates over the filenames and downloads each file using the 'retrbinary' command. Then, it selects
    the most recent version of the data for each available observation date.

    Note:
        - The 'secretcodes' module must contain the 'swot_username' and 'swot_password' variables with
          the appropriate FTP credentials.
    """

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

    os.makedirs(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+dataset_swot_n, exist_ok = True)
    os.chdir(currdir+'/input_'+name_experiment+'/'+today.strftime('%Y%m%d')+'/'+dataset_swot_n)

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