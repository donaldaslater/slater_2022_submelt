## Required code that is not included here

1. Topotoolbox (for flow routing). I used version 2.2 available from https://topotoolbox.wordpress.com
2. Additional colourbars available from https://www.mathworks.com/matlabcentral/fileexchange/28943-color-palette-tables-cpt-for-matlab

## Required datasets that are not included here

1. BedMachine Greenland, available from e.g. https://nsidc.org/data/idbmg4/versions/4. I used the version dated 2020-07-16
2. RACMO runoff data, available on request from Brice Noel at B.P.Y.Noel@uu.nl
3. Basal melt rate dataset from Karlsson 2021 as described and available at https://www.nature.com/articles/s41467-021-23739-z
4. Ocean reanalysis products ASTE, CHORE, EN4 and ORAS5. Instructions for downloading data for each of these are provided in /data/ASTE etc  

## Workflow for provided code

Note that the paths to the directories containing the above datasets will need to be changed to where you have put the datasets.

1. Run `data/get_drainagebasins.m`, which does the hydrological flow routing to define hydrological basins for the glaciers. This script also delineates the glacier calving fronts, which is a part automated/part manual process and may need to be updated if using a different version of BedMachine. The list of 243 major glaciers is taken from Morlighem2017_majorglaciers.xlsx, which was downloaded from the supporting information of the Morlighem BedMachine paper described above. I have manually added the BedMachine coordinates of each glacier to the spreadsheet. Ultimately this script creates drainagebasins.mat - a mask file where each pixel in a drainage basin is assigned the number of that glacier, and glaciers.mat - a matlab structure containing initial information about each of the glaciers.

2. Run `data/get_runoff.m`, which sums RACMO runoff over the hydrological drainage basin for each glacier. This adds the monthly runoff data to the glaciers.mat structure.

3. Run `data/get_basalmelting.m`, which sums the Karlsson basal melt rate dataset over each glacier's hydrological basin to get a subglacial discharge from basal melting, which is added to the glaciers.mat structure.

4. Run `data/get_oceanpoints.m`, which finds, for each glacier, the closest, deepest ocean product grid point. The (x,y,z) and linear index of this point is stored in the glaciers.mat structure and in the oceanindices.mat file. This script also produces plots of this procedure for sense checking.

5. Run the scripts `data/get_ocean_ASTE.py` or `data/get_ocean_CHORE.py` or `data/get_ocean_EN4.py` or `data/get_ocean_ORAS5.py`. These take the ocean index for each glacier from oceanindices.mat and use this to create a .mat file for each product (e.g. ocean_ASTE_glaciers.mat) that contains only the ocean profiles relevant to the glaciers. That is, ocean_ASTE_glaciers.mat contains only the ASTE grid points relevant to the glaciers, rather than all grid points. These scripts are written in python rather than matlab because, due to the volume of data, I had to run them on the server rather than my laptop, and the server I used does not have matlab.

6. Run `data/finalise_ocean.m`, which reads the .mat ocean files (e.g. ocean_ASTE_glaciers.mat) to associate the ocean profiles to the glacier in question. It also cuts off the ocean profiles for each glacier at the sill depth and interpolates to the get the ocean temperature and salinity time series at the grounding line of the glaciers. Saves in the glaciers.mat structure.

7. Run `data/finalise_glaciers.m`, which converts ocean temperature and salinity to ocean thermal forcing, combines time series from the different ocean products into one average time series and calculates submarine melt rate. Lastly, it filters glaciers to a final dataset used in the paper by removing glaciers that (i) have a grounding line depth >-50 m (which are either very small glaciers or have unreliable bed data), (ii) have a sill depth >-50 m (suggesting that the bathymetry is unreliable) and (iii) have a mean annual subglacial discharge <2.5 m3/s (suggesting these are very small glaciers or the hydrological basin is unreliable). These filters ultimately remove 120 glaciers, leaving a dataset of 123.

8. Run `model/glacier_response.m', which runs the simple glacier model forced by the submarine melt rate time series as described in the paper. It creates modeloutput.mat, which contains the simulation output required to create the plots.

9. Run 'plotting/makeplots.m' to create the plots for the main paper, or 'plotting/makeplots_extendeddata.m' to make the extended data figures.
