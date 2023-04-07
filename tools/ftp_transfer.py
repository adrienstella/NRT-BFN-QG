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
    ftp.cwd(month)

    # Set the name of the file to download
    filenames = ftp.nlst()
    for filename in filenames:
        print('Retreiving data for '+filename)

        # Download the file
        ftp.retrbinary("RETR "+filename, open(filename, 'wb').write)