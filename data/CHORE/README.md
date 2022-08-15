CHORE_RL ocean reanalysis is described in

Yang, C., Masina, S., and Storto, A.: Historical ocean reanalyses (1900-2010) using different data assimilation strategies, Quarterly Journal of the Royal Meteorological Society, pp. 479-493, https://doi.org/10.1002/qj.2936, 2017.

Downloading the data is a bit awkward: need to do  
`ftp downloads.cmcc.bo.it`  
with the username and password from http://c-glors.cmcc.it/index/index-7.html?sec=7 then  
`cd /p_cglors/CHOR/CHORE_RL`  
`prompt`  
`mget CHORE_RL_1979.nc.gz`  
`mget CHORE_RL_198*.nc.gz`  
`mget CHORE_RL_199*.nc.gz`  
`mget CHORE_RL_2*.nc.gz`  
and to get the grid file and landmask do  
`cd /p_cglors/CHOR/TOOLS`  
`get GRID.nc`  
then unpack the files using  
`gunzip *.gz`  
The unpacked files for 1979-2010 take up ~84GB.

Running the `get_oceanpoints.m` script requires only the GRID.nc file.
