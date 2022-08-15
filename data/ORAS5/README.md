ORAS5 is an ocean reanalysis run by ECMWF and described here: https://www.ecmwf.int/en/research/climate-reanalysis/ocean-reanalysis

Download data files using  
`wget -i ORAS5_downloadlist.txt`  
Unpack files using  
`cat *.tar.gz | tar zxvf - -i`  
Once unpacked, the ORAS5 monthly files for 1979-2018 occupy ~560GB.

Running the `get_oceanpoints.m` script just requires one of the output files, e.g. votemper_ORAS5_1m_201812_grid_T_02.nc
