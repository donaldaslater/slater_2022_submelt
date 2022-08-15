from scipy.io import netcdf
from scipy.io import savemat
from scipy.io import loadmat
import numpy as np
import os
import glob
import coords

# indices corresponding to glaciers
inds = loadmat("oceanindices.mat")
morlighem_n = inds['morlighem_n'].squeeze()
ind_ORAS5 = inds['ind_ORAS5'].squeeze()

# input files
filelistT = sorted(glob.glob(os.path.join('/scratch/other/das26/ORAS5/', 'votemper*.nc')))
filelistS = sorted(glob.glob(os.path.join('/scratch/other/das26/ORAS5/', 'vosaline*.nc')))
#filelistT = filelistT[0:24]
#filelistS = filelistS[0:24]

# loop over files to get time vector
t = np.empty(len(filelistT))
for i in range(len(filelistT)):
	datestr = filelistT[i].split("_")[3]
	t[i] = float(datestr[0:4]) + float(datestr[4:6])/12 - 1.0/24

# read depth
file2read = netcdf.NetCDFFile(filelistT[0],'r',version=1)
z = -file2read.variables['deptht'].data.astype('float64')

# initialise output arrays
potentialT = np.empty((len(ind_ORAS5),len(z),len(t)))
practicalS = np.empty((len(ind_ORAS5),len(z),len(t)))

# loop over netcdf files
for i in range(len(filelistT)):

	# read T/S
	file2readT = netcdf.NetCDFFile(filelistT[i],'r',version=1)
	Ti = file2readT.variables['votemper'].data.astype('float64').squeeze()
	file2readS = netcdf.NetCDFFile(filelistS[i],'r',version=1)
	Si = file2readS.variables['vosaline'].data.astype('float64').squeeze()
	Ti = np.reshape(Ti,(Ti.shape[0],Ti.shape[1]*Ti.shape[2]))
	Si = np.reshape(Si,(Si.shape[0],Si.shape[1]*Si.shape[2]))

	# set no data values to NaN
	nodata_inds = Ti>100
	Ti[nodata_inds] = float("nan")
	Si[nodata_inds] = float("nan")

	# assign to glaciers
	for j in range(len(ind_ORAS5)):
		if np.isnan(ind_ORAS5[j]):
			potentialT[j,:,i] = float("nan")*Ti[:,0]
			practicalS[j,:,i] = float("nan")*Si[:,0]
		else:
			potentialT[j,:,i] = Ti[:,ind_ORAS5[j]-1] # -1 to account of python starting at 0 versus matlab at 1
			practicalS[j,:,i] = Si[:,ind_ORAS5[j]-1]

# save
savemat("ocean_ORAS5_glaciers.mat",{"z":z,"potentialT":potentialT,"practicalS":practicalS,"t":t,"morlighem_n":morlighem_n},do_compression=True,oned_as='row')
