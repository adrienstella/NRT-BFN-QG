# NRT-BFN-QG

A program which produces near-real-time dynamic maps of ocean surface topography by assimilation of altimetric data into a quasi-geostrophic model.
See example/ for an example notebook with figures.

### Main steps of the algorithm: 
1) Download alongtrack (L3) altimetric data from a constellation of 7 nadirs on the CMEMS FTP server + SWOT nadir from AVISO
2) Download gridded (L4) DUACS ADT product on the CMEMS FTP server (For boundary conditions and mask)
3) Run the Back-and-Forth Nudging algorithm from MASSH (https://github.com/leguillf/MASSH) with the associated configuration, assimilating the data with the quasi-geostrophic model for the chosen spatial and temporal windows.
4) Process the outputs : cut off spin-up period, make daily averages, compute additional variables
5) Share created maps : copy .nc files to remote server (here, Ifremer)

Variables inclued : Sea Surface Height, Geostrophic (horizontal) velocities, Normalized relative vorticity.

### Parameters: 
- *lon_min*, *lon_max*, *lat_min*, *lat_max*: spatial domain
- *today*: final date for assimilation (default = current computer day and time)
- *numdays*: number of days of assimilation before final date

### Installation
:computer: _**How to get started ?**_

Clone this repo: 
```
git clone https://github.com/adrienstella/NRT-BFN-QG.git
```
Clone the MASSH repo:
```
git clone https://github.com/leguillf/MASSH.git
```

Adapt the different paths to match yours if needed
Enter your usernames and passwords for the FTP servers you want to connect to (or create your own secretcodes.py file).

You're now good to go ! Just run NRT_BFN_main.py to start mapping.
