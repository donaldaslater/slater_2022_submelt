## Required code that is not included here

## Required datasets that are not included here

1. BedMachine Greenland.
2. 


## Workflow

1. Download BedMachine Greenland topography dataset. The dataset is described in Morlighem et al, 2017, GRL, and can be downloaded by following links from here (https://sites.uci.edu/morlighem/dataproducts/bedmachine-greenland/) but I used a newer version (2020-07-16) obtained directly from Mathieu.

2. Run `get_drainagebasins.m`, which does the hydrological flow routing to define hydrological basins for the glaciers. This script also delineates the glacier calving fronts, which is a part automated/part manual process and may need to be updated if using a different version of BedMachine. The list of 243 major glaciers is taken from Morlighem2017_majorglaciers.xlsx, which was downloaded from the supporting information of the Morlighem BedMachine paper described above. I have manually added the BedMachine coordinates of each glacier to the spreadsheet. Ultimately this script creates drainagebasins.mat - a mask file where each pixel in a drainage basin is assigned the number of that glacier, and glaciers.mat - a matlab structure containing initial information about each of the glaciers.

3. Download RACMO data following, for example, instructions in the /RACMO folder.

4. Run `get_runoff.m`, which sums RACMO runoff over the hydrological drainage basin for each glacier. This adds the monthly runoff data to the glaciers.mat structure.

5. Download ice velocity data - the estimate of basal melting requires an ice velocity field for the ice sheet. I obtained this from https://catalogue.ceda.ac.uk/uuid/eaed9fba86c44e9c854dfbdec9d16b99, which is only a snapshot and a possible future improvement would be to make this time variable.

6. Run `get_basalmelting.m`, which does a simple calculation of basal melt rate based on a basal traction and ice velocity, and sums this over the glacier's hydrological basin to get a subglacial discharge from basal melting, which is added to the glaciers.mat structure.

7. Download enough of the ocean files to be able to load the grid and land mask for each ocean product considered. See individual ocean product directories (/ASTE, /CHORE, /EN4, /ORAS5) for details.

8. Run `get_oceanpoints.m`, which finds, for each glacier, the closest, deepest ocean product grid point. The (x,y,z) and linear index of this point is stored in the glaciers.mat structure and in the oceanindices.mat file. This script also produces plots of this procedure for sense checking.

9. Download full ocean product files following the instructions in the directories. This involves a large volume of data and so I did this on a server rather than locally.

10. Run the scripts `get_ocean_ASTE.py` or `get_ocean_CHORE.py` or `get_ocean_EN4.py` or `get_ocean_ORAS5.py`. These take the ocean index for each glacier from oceanindices.mat and use this to create a .mat file for each product (e.g. ocean_ASTE_glaciers.mat) that contains only the ocean profiles relevant to the glaciers. That is, ocean_ASTE_glaciers.mat contains only the ASTE grid points relevant to the glaciers, rather than all grid points. These scripts are written in python rather than matlab because, due to the volume of data, I had to run them on the server rather than my laptop, and the server does not have matlab.

11. Run `finalise_ocean.m`, which reads the .mat ocean files (e.g. ocean_ASTE_glaciers.mat) to associate the ocean profiles to the glacier in question. It also cuts off the ocean profiles for each glacier at the sill depth and interpolates to the get the ocean temperature and salinity time series at the grounding line of the glaciers. Saves in the glaciers.mat structure.

12. Run `finalise_glaciers.m`, which converts ocean temperature and salinity to ocean thermal forcing, combines time series from the different ocean products into one average time series and calculates submarine melt rate. Lastly, it filters glaciers to a final dataset used in the paper by removing glaciers that (i) have a grounding line depth >-50 m (which are either very small glaciers or have unreliable bed data), (ii) have a sill depth >-50 m (suggesting that the bathymetry is unreliable) and (iii) have a mean annual subglacial discharge <2.5 m3/s (suggesting these are very small glaciers or the hydrological basin is unreliable). These filters ultimately remove 120 glaciers, leaving a dataset of 123.
