from scipy.io import netcdf
from scipy.io import savemat
from scipy.io import loadmat
import numpy as np
import os
import glob
import coords

# input files
filelistT = sorted(glob.glob(os.path.join('/scratch/other/das26/ASTE/','THETA*.nc')))
filelistS = sorted(glob.glob(os.path.join('/scratch/other/das26/ASTE/','SALT*.nc')))

# read time and depth vectors and get lat/lon shape
file2read = netcdf.NetCDFFile(filelistT[0],'r',version=1)
z = -file2read.variables['dep'].data.astype('float64').squeeze()
lon = file2read.variables['lon'].data.astype('float64').squeeze()

# create manual t vector for ASTE
t = np.empty(12*16)
for i in range(len(t)):
	t[i] = 2002 + 1.0/24 + i*1.0/12

# initialise data arrays
T = np.empty((len(filelistT),len(t),len(z),lon.shape[0],lon.shape[1]))
S = np.empty((len(filelistS),len(t),len(z),lon.shape[0],lon.shape[1]))
lat = np.empty((len(filelistT),lon.shape[0],lon.shape[1]))
land = np.empty((len(filelistT),len(t),len(z),lon.shape[0],lon.shape[1]))
lon = np.empty((len(filelistT),lon.shape[0],lon.shape[1]))

# load data
for i in range(len(filelistT)):
	file2read = netcdf.NetCDFFile(filelistT[i],'r',version=1)
	T[i,:,:,:,:] = file2read.variables['THETA'].data.astype('float64').squeeze()
	lon[i,:,:] = file2read.variables['lon'].data.astype('float64').squeeze()
	lat[i,:,:] = file2read.variables['lat'].data.astype('float64').squeeze()
	file2read = netcdf.NetCDFFile(filelistS[i],'r',version=1)
	S[i,:,:,:,:] = file2read.variables['SALT'].data.astype('float64').squeeze()
	for j in range(len(t)):
		land[i,j,:,:,:] = file2read.variables['land'].data.astype('float64').squeeze()

# set land values to NaN
land_inds = land==0
T[land_inds] = float("nan")
S[land_inds] = float("nan")

# put into required shape
T = np.transpose(T,(0,4,3,2,1))
T = np.reshape(T,(T.shape[0]*T.shape[1]*T.shape[2],T.shape[3],T.shape[4]),order='F')
S = np.transpose(S,(0,4,3,2,1))
S = np.reshape(S,(S.shape[0]*S.shape[1]*S.shape[2],S.shape[3],S.shape[4]),order='F')

# indices corresponding to glaciers
inds = loadmat("oceanindices.mat")
morlighem_n = inds['morlighem_n'].squeeze()
ind_ASTE = inds['ind_ASTE'].squeeze()

# initialise output arrays
potentialT = np.empty((len(ind_ASTE),T.shape[1],T.shape[2]))
practicalS = np.empty((len(ind_ASTE),S.shape[1],S.shape[2]))

# assign to glaciers
for i in range(len(ind_ASTE)):
	if np.isnan(ind_ASTE[i]):
		potentialT[i,:,:] = float("nan")*T[0,:,:]
		practicalS[i,:,:] = float("nan")*S[0,:,:]
	else:
		potentialT[i,:,:] = T[ind_ASTE[i]-1,:,:] # -1 to account of python starting at 0 versus matlab at 1
		practicalS[i,:,:] = S[ind_ASTE[i]-1,:,:]

# save
savemat("ocean_ASTE_glaciers.mat",{"z":z,"potentialT":potentialT,"practicalS":practicalS,"t":t,"morlighem_n":morlighem_n},do_compression=True,oned_as='row')

