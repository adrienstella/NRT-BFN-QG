# NRT-BFN-QG

A program which produces near-real-time dynamic maps of ocean surface topography by assimilation of altimetric data into a quasi-geostrophic model.

### Main steps of the algorithm : 
1) Download alongtrack (L3) altimetric data from a constellation of 7 nadirs on the CMEMS FTP server
2) Download gridded (L4) DUACS ADT product on the CMEMS FTP server (For boundary conditions and mask)
3) Run the Back-and-Forth Nudging algorithm from MASSH (https://github.com/leguillf/MASSH) with the associated configuration, assimilating the data with the quasi-geostrophic model for the chosen spatial and temporal windows.
4) Process the outputs : cut off spin-up period, make daily averages, compute additional variables
5) Share created maps : copy .nc files to remote server (here, Ifremer)

Variables inclued : Sea Surface Height, Geostrophic (horizontal) velocities, Normalized relative vorticity.

** add more info on how to start
e.g. tell people to do a git clone from MASSH and remove the MASSH folder from here.
+ make a notebook (quickstart/example) that shows how everything runs with the outputs etc. with explanations
