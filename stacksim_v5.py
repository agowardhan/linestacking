import sys 
import time
import numpy as np
from astropy.io import fits
import scipy as sc  
import matplotlib.pyplot as plt
from scipy import spatial 
from pylab import  *	
from math import e 

c = 299792358.
hz2ghz = 10.**9

ion()



def uvsim(fits_input, file_out = 'uvsim_out.fits',  width=50, ln = [200, 500, 800, 1500], mod = 0.2, centre = 1150 ):

	
	'''

	Given a uvfits file, with a spectral line of width = <width>, centred on <centre>, this code insets the spectral linear at other channels, listed in <ln>. The amplitude is modulated by <mod>

	Inputs : 

	fits_input: string 
		enter input uvfits file	
	
	file_out: string
		output fits file name, default 'uvsim_out.fits'
	
	width: int 
		width of spectral line in channels. 
	
	ln: array([int]) 
		list of channels where spectral line is to be inserted 
	
	mod: float 
		amplitude scale from original line
	
	centre: int 
		channel where spectral line is centred 

	'''
	# channels 36 to 60, from 1215, to 1099 in the other one. 
	# data_cube, header_data_cube = fits.getdata(fits_input, 0, header=True)
	# file = fits.open('test.fits', mode='append')
	
	file = fits.open(fits_input)						# open uvfits file
	header1 = file[0].header  							# extract header 
	uvtable = file[0].data 								# extract data 
	vis = uvtable['DATA']								# extract visibilities from data 
	n_chan = np.shape(file[0])[3] 						# assuming shape of uvfits file to be something like (1, 1, 1, 100, 1, 3, 0)

	### beware - the shape may be different for different uvfits files 
	
	n_vis = np.shape(vis)[0]							# number of visibilities
	vis0 = vis 

# line width in channels - make sure its even. check what happens when its odd ?
	
	for chan in ln: 
		
		if ((chan < 100 + centre) and (chan > centre - 100)): 
			print" Beware : adding spectral line within 200 channels of the original.."

		for i in range(0,n_vis):									# for each visibility, stacking the spectral lines. 

			tmp1 = vis[i][0][0][0][chan-width/2:chan+width/2]		# extracted original weights - keeping them the same. 
			wts = tmp1[:,:,2]

			tmp2 = vis[i][0][0][0][centre-width/2:centre +width/2] 	# Extracting spectral line 
			re = tmp2[:,:,0]
			im = tmp2[:,:,1]

			tmp3 = np.zeros((width,1,3))
			tmp3[:,:,0] = re
			tmp3[:,:,1] = im
			tmp3[:,:,2] = wts 
			#tmp3[:,:,2] = np.zeros((np.shape(re))) 

			vis[i][0][0][0][chan-width/2:chan+width/2]+= np.multiply(tmp3,mod)


	uu = uvtable['UU']				# copying over all the data to create the new fits file 
	vv = uvtable['VV']				
	ww = uvtable['WW']
	base = uvtable['BASELINE']
	date1 = uvtable['DATE']
	date2 = uvtable['_DATE']

	# uu_given = (uu_orig + pzero)* pscal 
	# unscaling and unshifting the given u,v,w visibilities. 

	pars = [uu,vv,ww,base2,date21,date22]

 	pscal = header1['PSCAL*']     
 	pzero = header1['PZERO*']     
	
	for i in range(0,5): 											# shifting and scaling parameters 
		pars[i] = [x - pzero[i] for x in pars[i]]
		pars[i]/= pscal[i]

	pnames=['UU','VV','WW','BASELINE','DATE','DATE']
	pdata=pars
	
	tmp4= fits.GroupData(vis,bitpix=-32, pardata=pdata, parnames=pnames)	# making data structure 
	hdu = fits.GroupsHDU(tmp4, header1)										# attaching header 
	tmp5 = fits.HDUList([hdu, file[1]])										# creating HDU 
 	tmp5.writeto('temp_sim.fits')  # the GroupdHDU must be the first one. 			# writing out file. 

 	tmp6 = fits.open('temp_sim.fits')
 	tmp6[0].header['BSCALE'] = header1['BSCALE']
 	tmp6[0].header['BZERO'] = header1['BZERO']
 	tmp6.writeto(file_out, clobber=clob)
