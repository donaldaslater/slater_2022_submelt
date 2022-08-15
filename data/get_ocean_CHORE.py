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
ind_CHORE = inds['ind_CHORE'].squeeze()

# input files
gridfile = sorted(glob.glob(os.path.join('/scratch/other/das26/CHORE/','GRID*.nc')))
filelist = sorted(glob.glob(os.path.join('/scratch/other/das26/CHORE/','CHORE*.nc')))
#filelist = filelist[0:5]

# load grid and mask
file2read = netcdf.NetCDFFile(gridfile[0],'r',version=1)
z = -file2read.variables['dep'].data.astype('float64')
landmask0 = file2read.variables['tmsk'].data.astype('float64')
landmask = np.empty((12,landmask0.shape[0],landmask0.shape[1],landmask0.shape[2]))
for i in range(12):
	landmask[i,:] = landmask0

# month and time vector
m = np.empty(12)
for i in range(len(m)):
	m[i] = 1.0/24 + i*1.0/12
t = np.empty(12*len(filelist))
for i in range(len(filelist)):
	yr = float(filelist[i].split("_")[2].split(".")[0])
	t[12*i:12*(i+1)] = yr+m

# initialise output arrays
potentialT = np.empty((len(ind_CHORE),len(z),len(t)))
practicalS = np.empty((len(ind_CHORE),len(z),len(t)))

# loop over netcdf files
for i in range(len(filelist)):

	# read T/S
	file2read = netcdf.NetCDFFile(filelist[i],'r',version=1)
	Ti = file2read.variables['votemper'].data.astype('float64').squeeze()
	Ti = Ti*file2read.variables['votemper'].scale_factor+file2read.variables['votemper'].add_offset
	Si = file2read.variables['vosaline'].data.astype('float64').squeeze()
	Si = Si*file2read.variables['vosaline'].scale_factor+file2read.variables['vosaline'].add_offset
	# mask for land
	land_inds = landmask==0
	Ti[land_inds] = float("nan")
	Si[land_inds] = float("nan")
	# put into required shape
	Ti = np.transpose(Ti,(3,2,1,0))
	Si = np.transpose(Si,(3,2,1,0))
	Ti = np.reshape(Ti,(Ti.shape[0]*Ti.shape[1],Ti.shape[2],Ti.shape[3]),order='F')
	Si = np.reshape(Si,(Si.shape[0]*Si.shape[1],Si.shape[2],Si.shape[3]),order='F')

	# assign to glaciers
	for j in range(len(ind_CHORE)):
		if np.isnan(ind_CHORE[j]):
			potentialT[j,:,12*i:12*(i+1)] = float("nan")*Ti[0,:,:]
			practicalS[j,:,12*i:12*(i+1)] = float("nan")*Si[0,:,:]
		else:
			potentialT[j,:,12*i:12*(i+1)] = Ti[ind_CHORE[j]-1,:,:] # -1 to account of python starting at 0 versus matlab at 1
			practicalS[j,:,12*i:12*(i+1)] = Si[ind_CHORE[j]-1,:,:]

# save
savemat("ocean_CHORE_glaciers.mat",{"z":z,"potentialT":potentialT,"practicalS":practicalS,"t":t,"morlighem_n":morlighem_n},do_compression=True,oned_as='row')
