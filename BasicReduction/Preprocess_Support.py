# This program provides support functions for Preprocess module.

# Basic library imports.

import numpy as np
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import stsdas

import pyfits as pft
import os,sys,glob,string,time
import math as mth
import pickle

irafhome = '/home/kaustubh/iraf/'


# Basic input functions.
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

# Function for running lacosmic, a cosmic ray removal algorithm.
def lacosmic(name,irafhome):
	import time
	iraf.task(lacos_spec=irafhome+'lacos_spec.cl')
	outname = 'la'+name
	pl = 'mask'+name
	iraf.lacos_spec(input=name,output=outname,outmask=pl,gain=1.,readn=2.89,
				xorder=9,yorder=0,sigclip=4.5,sigfrac=0.5,objlim=1.,niter=7,verbose='yes',mode='al')
	old = time.time()
	os.system('ds9 %s -zscale %s -zscale -blink &' % (pl,name))
	os.system('ds9 %s -zscale &' % outname)
	return

