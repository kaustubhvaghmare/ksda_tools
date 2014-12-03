# Support for file for Flat.py

# Import libraries.
import numpy as np
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import stsdas

import pyfits as pft
import os,sys,glob,string,time
import math as mth
import pickle

# Import supporting functions
from Flat_Support import *

# Define the IRAF home directory.
irafhome = '/home/kaustubh/iraf/'

# Basic input / output functions.
def input_val(script,lim1,lim2):
	while True:
	   try:
	       value = int(raw_input('%s : ' % (script)))
	   except ValueError: # just catch the exceptions you know!
	       print 'That\'s not a number!'
	   else:
	       if lim1 <= value < lim2: 
	           break
	       else:
	           print 'Out of range. Try again'
	return value

#------------------------------------------------------------------------------------------------------------------------------------

def input_str(script):
	while True:
		userInput = str(raw_input('%s' % (script)))
		if len(userInput) == 1:
			if userInput in string.letters:
				if userInput.lower() == 'y' or userInput.lower() == 'n':		
					break
				print 'Please enter only "y" or "n"'
			else:
				print 'Please enter only letters!'
		elif len(userInput) == 0:
			print 'Please enter at least 1 character!'
		elif len(userInput) >1 and userInput.isalpha():
			print 'Please enter only 1 character!'
		else:
			print 'Please enter only letters and no numbers'
	return userInput.lower()

# The flat fielding function

def flat(flist):
	f = open('flatlist','w')
	for i in range(0,len(flist)):
		f.write('%s\n' % (flist[i]))
	f.close()
	f_list = '@flatlist'
	combine_flat = 'com_flat.fits'

	print "Combining images."
	iraf.images.immatch.imcombine(input=f_list,output=combine_flat,headers='',bpmasks='',rejmask='',nrejmas='',expmask='',
					sigmas ='',imcmb='$I',logfile='STDOUT',combine='median',reject ='none',project='no',
					outtype='real',outlimi='',offsets='none',masktyp='none',maskval=0,blank=0.,scale='none',
					zero='none',weight='none',statsec='',expname='',lthresh='INDEF',hthresh='INDEF',nlow=1,nhigh=1,
					nkeep=1,mclip='yes',lsigma=3.,hsigma =3.,rdnoise=0.,gain=1.,snoise =0.,sigscal=0.1,
					pclip=-0.5,grow=0.,mode='ql')
	
	print "Filling gaps in CCD. Needed to avoid artifical gradients while performing illumination correction."	
	ccdgap(combine_flat) # This produces a file called ccom_flat.fits. This needs 
	combine_flat = "c" + combine_flat
	illum_flat = 'il'+combine_flat

	iraf.noao.imred.ccdred.ccdproc.noproc='no'
	iraf.noao.imred.ccdred.ccdproc.fixpix='no'
	iraf.noao.imred.ccdred.ccdproc.oversca='no'
	iraf.noao.imred.ccdred.ccdproc.trim='no'
	iraf.noao.imred.ccdred.ccdproc.zerocor='no'
	iraf.noao.imred.ccdred.ccdproc.darkcor='no'
	iraf.noao.imred.ccdred.ccdproc.flatcor='no'
	iraf.noao.imred.ccdred.ccdproc.illumco='yes'
	iraf.noao.imred.ccdred.ccdproc.fringec='no'
	iraf.noao.imred.ccdred.ccdproc.readcor='no'
	iraf.noao.imred.ccdred.ccdproc.scancor='no'

	"Performing illumination correction on Flat."
	iraf.noao.imred.ccdred.mkillumflat(input=combine_flat,output=illum_flat,ccdtype='',xboxmin=3.,xboxmax=5,yboxmin=3.,
					yboxmax=5,clip='yes',lowsigm=2.5,highsig=2.5,divbyze=1.,ccdproc='',mode='ql')
	cillum_flat = illum_flat
#	os.system('ds9 %s &' % (cillum_flat))
	fimg = pft.open(cillum_flat)
	prihdr = fimg[0].header
	scidata = fimg[0].data
	n1 = prihdr['NAXIS1']
	n2 = prihdr['NAXIS2']

	print "Normalizing Flats using Median of Entire Frame."
	
	median_flux = np.median(scidata[0:n2,0:n1])
	scidata1 = scidata/median_flux
	new_name = 'master_flat.fits'
	pft.writeto(new_name,data=scidata1,header=prihdr,clobber=True)
	for i in range(0,len(flist)):
		os.system('mv %s history/' % (flist[i]))
	os.system('mv %s history/' % (combine_flat))
	os.system('mv %s history/' % (illum_flat))
	os.system('mv %s history/' % (cillum_flat))
	return

# Performs division.
def divide(name1,name2,switch):
	ffts1 = pft.open(name1)
	prihdr1 = ffts1[0].header
	scidata1 = ffts1[0].data
	ffts2 = pft.open(name2)
	prihdr2 = ffts2[0].header
	scidata2 = ffts2[0].data
	x,y  = np.where(scidata2==0)
	for i in range(0,len(x)):
		scidata2[x[i],y[i]] = 1.
	scidata3 = abs(scidata1/scidata2)
	if switch == 1:
		new = 'd'+name1
	else:
		new = 'fl'+name1

	ffts3 = pft.writeto(new,scidata3,header=prihdr1,clobber=True)
	ffts1.close()
	ffts2.close()
	#os.system('mv %s history/' % (name1))
	#os.system('mv %s history/' % (name2))
	return

# The ccdgap interpolation function

def ccdgap(name):
	fimg = pft.open(name)
	prihdr = fimg[0].header
	scidata = fimg[0].data

	n1 = prihdr['NAXIS1']
	n2 = prihdr['NAXIS2']

	#below are the 4 coordinates of edge of the ccd gap
	a = ccd_locate(scidata)[0]-4;b = ccd_locate(scidata)[1]+2;c = ccd_locate(scidata)[2]-2;d = ccd_locate(scidata)[3]+3
	e = n2 #ccd height

	gap1_part1 = (scidata[:,a-6:a-1].sum(axis=1))/5.0
	gap1_part2 = (scidata[:,b+1:b+6].sum(axis=1))/5.0
	gap2_part1 = (scidata[:,c-6:c-1].sum(axis=1))/5.0
	gap2_part2 = (scidata[:,d+1:d+6].sum(axis=1))/5.0
	
	grad1 = (gap1_part2-gap1_part1)/((b-a)+5.0)
	grad2 = (gap2_part2-gap2_part1)/((d-c)+5.0)

	for i in range(a,b):
		scidata[:,i] = grad1*((i-a)+2)+gap1_part1

	for i in range(c,d):
		scidata[:,i] = grad2*((i-c)+2)+gap2_part1

	namec = "c"+name
	pft.writeto(namec,data=scidata,header=prihdr,clobber=True)
	fimg.close()
	os.system('mv %s history/' % (name))
	return

#function to help locate the ccd gap automatically
def ccd_locate(data):
	sum_col = data.sum(axis=0) # summing over each row of the image
	first_gap = np.where(sum_col[0:len(sum_col)/2] == 0)[0]
	second_gap = (len(sum_col)/2)+np.where(sum_col[len(sum_col)/2:] == 0)[0]
	return first_gap[0],first_gap[-1],second_gap[0],second_gap[-1]


