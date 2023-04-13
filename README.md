# NRT-BFN-QG

A program which produces near-real-time dynamic maps of ocean surface topography by assimilation of altimetric data into a quasi-geostrophic model.

### Main steps of the algorithm: 
1) Download alongtrack (L3) altimetric data from a constellation of 7 nadirs on the CMEMS FTP server
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

Clone the data challenge repo: 
```
git clone https://github.com/adrienstella/NRT-BFN-QG.git
```
Enter the NRT-BFN-QG repo and clone the MASSH repo inside it:
```
git clone https://github.com/leguillf/MASSH.git
```

You're now good to go ! Just run NRT_BFN_main_2.0.py to start creating maps.




*Note to self: make a notebook (quickstart/example) that shows how everything runs with the outputs etc. with explanations*
