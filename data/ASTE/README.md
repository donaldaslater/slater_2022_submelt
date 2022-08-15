The Arctic State Estimate (ASTE) is described in Nguyen, A. T., Pillar, H., Ocaa, V., Bigdeli, A., Smith, T. A. et al. The Arctic Subpolar Gyre
sTate Estimate: Description and Assessment of a Data-Constrained, Dynamically Consistent Ocean-Sea Ice Estimate for 2002{2017. Journal of Advances in Modeling Earth Systems 13, e2020MS002398 (2021).

To visualise the tiles, I downloaded the sea surface height climatology for all tiles from https://web.corral.tacc.utexas.edu/OceanProjects/ASTE/Release1/nctiles_climatology/ETAN/ - I chose this variable simply because this is a small dataset that allows us to visualise the files.

Run `gridding.m` to visualise the tiles in relation to the coast of Greenland (coastline.mat). From this I identified that I needed tiles 5,11,12,14,15,27.

On the command line, run `wget -i ASTE_downloadlist.txt` to download the required T/S files.

To get the ocean points using `get_oceanpoints.m` requires a local copy of the T climatology files (i.e. T files without the time axis), which can be obtained from https://web.corral.tacc.utexas.edu/OceanProjects/ASTE/Release1/nctiles_climatology/THETA/
