The Arctic State Estimate (ASTE) is produced by An Nguyen and Patrick Heimach but has not yet been published. The data is arranged in tiles, so it is first necessary to decide which tiles are needed.

To visualise the tiles, I downloaded the sea surface height climatology for all tiles from https://web.corral.tacc.utexas.edu/OceanProjects/ASTE/Release1/nctiles_climatology/ETAN/ - I chose this variable simply because this is a small dataset that allows us to visualise the files.

Run `gridding.m` to visualise the tiles in relation to the coast of Greenland (coastline.mat). From this I identified that I needed tiles 5,11,12,14,15,27.

On the command line, run `wget -i ASTE_downloadlist.txt` to download the required T/S files.

To get the ocean points using `get_oceanpoints.m` requires a local copy of the T climatology files (i.e. T files without the time axis), which can be obtained from https://web.corral.tacc.utexas.edu/OceanProjects/ASTE/Release1/nctiles_climatology/THETA/
