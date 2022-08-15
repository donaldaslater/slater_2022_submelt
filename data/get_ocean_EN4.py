from scipy.io import netcdf
from scipy.io import savemat
from scipy.io import loadmat
import numpy as np
import os
import glob

# indices corresponding to glaciers
inds = loadmat("oceanindices.mat")
morlighem_n = inds['morlighem_n'].squeeze()
ind_EN4 = inds['ind_EN4'].squeeze()

# input files
filelist = sorted(glob.glob(os.path.join('/scratch/other/das26/EN4/', '*.nc')))
#filelist = filelist[0:2]

# loop over files to get time vector
t = np.empty(len(filelist))
for i in range(len(filelist)):
	datestr = filelist[i].split(".")[7]
	t[i] = float(datestr[0:4]) + float(datestr[4:6])/12 - 1.0/24

# read depth
file2read = netcdf.NetCDFFile(filelist[0],'r',version=2)
z = -file2read.variables['depth'].data.astype('float64')

# initialise output arrays
potentialT = np.empty((len(ind_EN4),len(z),len(t)))
practicalS = np.empty((len(ind_EN4),len(z),len(t)))

# loop over netcdf files
for i in range(len(filelist)):
	# read T/S
	file2read = netcdf.NetCDFFile(filelist[i],'r',version=2)
	Ti = file2read.variables['temperature'].data.astype('float64').squeeze()
	Si = file2read.variables['salinity'].data.astype('float64').squeeze()
	# account for scale and offset (would be automatic if using netcdf4)
	Ti = Ti*file2read.variables['temperature'].scale_factor+file2read.variables['temperature'].add_offset
	Si = Si*file2read.variables['salinity'].scale_factor+file2read.variables['salinity'].add_offset
	# restrict to Greenland
	Ti = np.reshape(Ti,(Ti.shape[0],Ti.shape[1]*Ti.shape[2]))
	Si = np.reshape(Si,(Si.shape[0],Si.shape[1]*Si.shape[2]))

	# convert to celsius
	Ti = Ti-273.15

        # set no data values to NaN
        nodata_inds = Ti<-100
        Ti[nodata_inds] = float("nan")
        Si[nodata_inds] = float("nan")

	# assign to glaciers
	for j in range(len(ind_EN4)):
		if np.isnan(ind_EN4[j]):
			potentialT[j,:,i] = float("nan")*Ti[:,0]
			practicalS[j,:,i] = float("nan")*Si[:,0]
		else:
			potentialT[j,:,i] = Ti[:,ind_EN4[j]-1] # -1 to account of python starting at 0 versus matlab at 1
			practicalS[j,:,i] = Si[:,ind_EN4[j]-1]

# save
savemat("ocean_EN4_glaciers.mat",{"z":z,"potentialT":potentialT,"practicalS":practicalS,"t":t,"morlighem_n":morlighem_n},do_compression=True,oned_as='row')

