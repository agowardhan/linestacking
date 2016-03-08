import sys 
import time
import numpy as np
from astropy.io import fits
import scipy as sc  
import matplotlib.pyplot as plt
from scipy import spatial 
from pylab import  *	
from math import e 

hz2ghz = 10.**9

ion()

def uvstack(fits_input, ln = [], freqs = [], width=1, file_out = 'uvs_out.fits', clob=False):

	'''

	Given a uvfits file, spectral line frequencies/channels and expected line width stacks spectral lines for each visibility

	Inputs: 

	fits_input: string 
		enter input uvfits file e.g. 'input.fits'
	
	file_out: string 
	 	output fits file name e.g. uvsim_out.fits
	
	freqs: float 
		line frequencies in GHz (default = [])
	
	width: int
	 	width of spectral line in channels. 
	
	ln: array[int] 
		list of channel number where spectral lines are centred 
	
	clob: bool 
		To overwrite existing file with same name as output  

	'''

	#
	# data_cube, header_data_cube  = fits.getdata(fits_input, 0, header=True)
	# file = fits.open('test.fits', mode='append')

	file = fits.open(fits_input)									# opening uvfits file
	nfits = np.shape(file)[0] 
	flag2 = 0 
	if (nfits > 1): 
		print'Not a UVFITS file:',nfits,'HDUs present' 
		print'Will append them to FITS file'
		flag2 = 2 
	header1  = file[0].header 										# extracting header for HDU 

	ln, freqs, flag = chkvals(header1, ln, freqs, width)			# function to check if frequencies/channel numbers are consistent with the FITS file 

	freq0 = header1['CRVAL4']										# extracting the central frequency
	n_chan = header1['NAXIS4']										# extracting central channel value 

	if(flag < 1): 
		return 

	uvtable = file[0].data  										# extracting uvtable  
	vis = uvtable['DATA']											# extracting visibiities 
	n_vis = np.shape(vis)[0]
	n_lines = np.shape(ln)[0]										# number of lines to stack

# line width in channels - make sure its even. 

	u_cent = uvtable['UU']											# U, V, W values given based on central frequency. 
	v_cent = uvtable['VV']
	w_cent = uvtable['WW']
	base = uvtable['BASELINE']
	date1 = uvtable['DATE']
	date2 = uvtable['_DATE']

	new_uvtable = np.zeros((n_vis*n_lines,))						# new empty table in shape of stacked visibilities 
	new_vis = np.zeros((n_vis*n_lines,1,1,1,width,1,3))             # this is FILE-SPECIFIC - beware.  

	# empty u,v,w,baseline, date arrays, 

	uu = []
	vv = [] 
	ww = [] 
	base2 = [] 
	date21 = [] 
	date22 = [] 

	for fq in freqs: 
		uu.extend(np.multiply(freq0/fq, u_cent))					# appending scaled u,v,w values to the empty visi file 
		vv.extend(np.multiply(freq0/fq, v_cent))
		ww.extend(np.multiply(freq0/fq, w_cent))
		base2.extend(base)
		date22.extend(date2)
		date21.extend(date1)
	
	for i in range(n_vis*n_lines):								# for all visibilities, copying over select channels 
		j = i/n_vis
		jj = i%n_vis  
		low, high = ln[j] - width/2, ln[j] + width/2 
		if (width%2!=0): 
			high = high + 1  										# fixing the shape for odd width
		if(flag2 > 0): 												# file format is not uv 
			tmp1 = vis[jj,:,:,:,low:high,:,:] 
		else: 
			tmp1 = vis[jj,:,:,low:high,:,:] 
		new_vis[i] = tmp1 				 # appends to the new visibility

	pars = [uu,vv,ww,base2,date21,date22]

 	pscal = header1['PSCAL*']     
 	pzero = header1['PZERO*']     
	
	for i in range(0,5): 											# shifting anc scaling parameters 
		pars[i] = np.array([x - pzero[i] for x in pars[i]])
		pars[i]/= pscal[i]

	# creating new GroupData

	pnames=['UU','VV','WW','BASELINE','DATE','DATE']   # note that 2 DATE names will result in one DATE and one _DATE columns. 
	pdata=pars
	tmp2= fits.GroupData(new_vis,bitpix=-32,  pardata=pdata, parnames=pnames) # bitpix -64 will not be understood by GILDAS. 
	hdu = fits.GroupsHDU(tmp2, header1)


	# writing data to fits file
	hdu_array = file
	hdu_array[0] = hdu
	tmp3 = fits.HDUList(hdu_array)
 	tmp3.writeto('temp.fits', clobber=True)  # the GroupdHDU must be the first one. 

 	# patch-job to include BSCALE/BZERO values which are not included in copying the header over. 

 	tmp4 = fits.open('temp.fits')
 	tmp4[0].header['BSCALE'] = header1['BSCALE']
 	tmp4[0].header['BZERO'] = header1['BZERO']
 	tmp4.writeto(file_out, clobber=clob)

# 	c = fits.open('uvs_out.fits')
# 	print c[0].data

##################
# IMAGE STACKING # 
##################

def imstack(fits_input, file_out = 'imgs_out.fits',  width=1, ln = [137,267,611,683, 842] ,clob=True, plt=0): 


	'''

Given an image fits file, spectral line frequencies/channels and expected line width, it stacks spectral lines on the image plane. 

	Inputs: 

	fits_input: string 
		enter input uvfits file e.g. 'input.fits'
	
	file_out: string 
	 	output fits file name e.g. 'uvsim_out.fits'
	
	freqs: float 
		line frequencies in GHz (defaule = [])
	
	width: int
	 	width of spectral line in channels. 
	
	ln: array[int] 
		list of channel number where spectral lines are centred 
	
	plt: int [default = 0]
		To plot stacked image or not

	'''
	flag = 2 
	img_cube, header  = fits.getdata(fits_input, 0, header=True)
	n_chan = np.shape(img_cube)[1]
	xl = np.shape(img_cube)[2]
 	yl = np.shape(img_cube)[3]

 	flag = chkvals_img(ln,width, n_chan)

 	if (flag < 1): 
		return 
	img_f = np.zeros((1,width,xl,yl))

	lns = np.shape(ln)[0]

	for l in ln:  
		img_f = img_f + img_cube[:,l-width/2.:l+width/2,:,:]
	img_f = img_f/lns
	if (plt > 0): 
		implt = plt.imshow(np.mean(img_f, axis=0))
		plt.colorbar()
		plt.show()

	header['NAXIS3'] = (width/2) + (width/2)
	fits.writeto(file_out, img_f, header, clobber=clob)

# function to extract spectra line information from uvfits fiel header

def chkvals(header, ln, freqs, width):


	'''

	Checks consistency between line channels, line frequencies, and values in header. 

	Inputs: 

	header: FITS header object  
		header from PrimaryHDU with visibility data 

	freqs: float 
		line frequencies in GHz (defaule = [])
	
	width: int
	 	width of spectral line in channels. 
	
	ln: array[int] 
		list of channel number where spectral lines are centred 
	
	'''

	flag = 2 
	freq0 = header['CRVAL4']
	n_chan = header['NAXIS4']
	ref =  header['CRPIX4']
	freq_res = header['CDELT4']
	freq_low = freq0 - ref*freq_res
	wd = freq_res*n_chan

	print '\n' , 'Central freq (GHz): ', freq0/hz2ghz, '\n' ,'Nchan:',  n_chan, '\n' ,'Ref. channel: ', ref, '\n' ,'Freq. res. (Hz):', freq_res, '\n','Lower freq(Hz): ',  freq_low, '\n' ,'Bandwidth (Hz): ', wd ,'\n' 
	
	for i in range(np.shape(ln)[0]-1):
		if (ln[i] > ln[i+1]):
			print('Channels unsorted. Getting out! Channels must be in order. Try again. ')
			flag = -1 

	for i in range(np.shape(ln)[0]-1): 
		if (ln[i] + width/2 > ln[i+1] - width/2):
			print'Warning - overlap between stacks',i,i+1

	if (ln==[] and freqs == []): 
		print '\n' ,'Channels/ frequencies of expected lines not provided - will use defaults'
		ln = [137,267,611,683,842]
		freqs=np.multiply(hz2ghz, [80.683, 81.190, 82.533, 82.815, 83.436])
		print'\n' ,'channels' , ln
		print'\n' ,'frequencies (in GHz)', freqs
		# will assume the HCN lines if none are given 

	elif (ln == []): 
		print"Channels of expected lines not provided - will calculate using frequencies"
		ln = [int(round((i - freq_low)*n_chan/wd)) for i in freqs] 
		print'Found channels', ln, 'from frequencies', freqs

	elif (freqs == []): 
		print"Frequencies of expected lines not provided - calculating channel numbers using frequencies"
		freqs = [((1.*i/n_chan)*wd + freq_low) for i in ln] 
		print'Calculated frequencies (Hz)', freqs, 'from channels' , ln

	if(np.max(ln) + width/2 > n_chan or np.min(ln) - width/2 < 0): 
		print'Invalid channel, width or frequency list given, check input'
		print(np.max(ln), n_chan, width)
		flag = -1 
	
	return ln, freqs, flag


def chkvals_img(ln, width, n_chan):


	'''

	Checks consistency between line channels and line width for image stacking.  

	Inputs: 
	
	width: int
	 	width of spectral line in channels. 
	
	ln: array[int] 
		list of channel number where spectral lines are centred 
	
	'''
	flag = 0 
	for i in range(np.shape(ln)[0]-1):
		if (ln[i] > ln[i+1]):
			print('Channels unsorted. Getting out! Channels must be in order. Try again. ')
			flag = -1 
			
	for i in range(np.shape(ln)[0]-1): 
		if (ln[i] + width/2 > ln[i+1] - width/2):
			print'Warning - overlap between stacks',i,i+1

	if(np.max(ln) + width/2 > n_chan or np.min(ln) - width/2 < 0): 
		print'Invalid channel, width or frequency list given, check input'
		print(np.max(ln), n_chan, width)
		flag = -1 
	
	return flag


def uvsim(fits_input, file_out = 'uvsim_out.fits',  width=50, ln = [200, 500, 800, 1500], mod = 0.2, centre = 1150, clob=True):

	
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

	pars = [uu,vv,ww,base,date1,date2]

 	pscal = header1['PSCAL*']     
 	pzero = header1['PZERO*']     
	
	for i in range(0,5): 											# shifting and scaling parameters 
		pars[i] = np.array([x - pzero[i] for x in pars[i]])
		pars[i]/= pscal[i]

	pnames=['UU','VV','WW','BASELINE','DATE','DATE']
	pdata=pars
	
	tmp4= fits.GroupData(vis,bitpix=-32, pardata=pdata, parnames=pnames)	# making data structure 
	hdu = fits.GroupsHDU(tmp4, header1)										# attaching header 
	tmp5 = fits.HDUList([hdu, file[1]])										# creating HDU 
 	tmp5.writeto('temp_sim.fits', clobber=True)  # the GroupdHDU must be the first one. 			# writing out file. 

 	tmp6 = fits.open('temp_sim.fits')
 	tmp6[0].header['BSCALE'] = header1['BSCALE']
 	tmp6[0].header['BZERO'] = header1['BZERO']
 	tmp6.writeto(file_out, clobber=clob)


def img_avg(fits_input, file_out = 'imgs_out.fits', width = -1,  ln = [-1, -1], clob = 0): 


	'''

Given an image fits file, spectral line frequencies/channels and expected line width, it stacks spectral lines on the image plane. 

	Inputs: 

	fits_input: string 
		enter input uvfits file e.g. 'input.fits'
	
	file_out: string 
	 	output fits file name e.g. 'uvsim_out.fits'
	
	freqs: float 
		line frequencies in GHz (defaule = [])
	
	width: int
	 	width of spectral line in channels. 
	
	ln: array[int] 
		list of channel number where spectral lines are centred 

	'''
	
	img_cube, header  = fits.getdata(fits_input, 0, header=True)
	nchan = np.shape(img_cube)[1]

	xl = np.shape(img_cube)[2]
 	yl = np.shape(img_cube)[3]

	if ((width < 0 or width > nchan) and ln[0]<0 or ln[1]< 0): 
		print 'Error - no valid range provided'
		return 
	elif (width < 0):
		nchanf = 1 
		img_f = np.mean(img_cube[:,ln[0]:ln[1],:,:],1) # averaging along axis 1 
		flag = 1 
	elif(ln == []): 
		nchanf = nchan/width + 1
		img_f = np.zeros((1,nchanf,xl,yl))
		for i in range(0, nchanf-1): 
			img_f[:,i,:,:] = np.mean(img_cube[:,i*width:(i+1)*width,:,:],1) # averaging along axis 1 
		# img_f[nchanf] = np.mean(img_cube[:,nchanf*width:nchan,:,:],1)
		flag = 2  		

	print(np.shape(img_f))
	
	if (plt > 0): 
		implt = plt.imshow(np.mean(img_f, axis=0))
		plt.colorbar()
		plt.show()

	header['NAXIS3'] = nchanf
	fits.writeto(file_out, img_f, header, clobber = clob)
	return img_f 

# function to extract spectra line information from uvfits fiel header