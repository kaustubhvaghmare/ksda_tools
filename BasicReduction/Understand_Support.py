
import numpy as np
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import stsdas

import pyfits as pft
import os,sys,glob,string,time
import math as mth
import pickle

from IRAFHome import *



# Support functions for Understand.py
def rename(fits_name):
	try:
		fimg = pft.open(fits_name)
		prihdr = fimg[0].header
		prihdr1 = fimg[1].header
		scidata = fimg[1].data

		n = prihdr1['NAXIS']
		n1 = prihdr1['NAXIS1']
		n2 = prihdr1['NAXIS2']
		object_name = prihdr['OBJECT']
		split_name = object_name.split()
		new_name = ''
		for i in range(0,len(split_name)):
			if i == 0: new_name = new_name+split_name[i]
			else: new_name = new_name+'_'+split_name[i]

		index = fits_name.split('.')[0][-3:]
		if len(new_name)>=11: new_name = new_name[:11]+'_'+index+'.fits'
		else: new_name = new_name+'_'+index+'.fits'

		telalt = prihdr['TELALT']
		airmass = 1.0/mth.sin(mth.radians(telalt))
		prihdr.update('airmass', airmass)

		prihdr.update('NAXIS', n);prihdr.update('NAXIS1', n1);prihdr.update('NAXIS2', n2)
		pft.writeto(new_name,data=scidata,header=prihdr,clobber=True)
		out = 1
	except:
		fimg = pft.open(fits_name)
		prihdr = fimg[0].header
		scidata = fimg[0].data

		object_name = prihdr['OBJECT']
		split_name = object_name.split()
		new_name = ''
		for i in range(0,len(split_name)):
			if i == 0: new_name = new_name+split_name[i]
			else: new_name = new_name+'_'+split_name[i]
	
		index = fits_name.split('.')[0][-3:]
		if len(new_name)>=11: new_name = new_name[:11]+'_'+index+'.fits'
		else: new_name = new_name+'_'+index+'.fits'

		telalt = prihdr['TELALT']
		airmass = 1.0/mth.sin(mth.radians(telalt))
		prihdr.update('airmass', airmass)

		pft.writeto(new_name,data=scidata,header=prihdr,clobber=True)
		out = 2
	return out

#function to read the header of fits files and output a few critical info needed in the program
def info_fits(fits_name):
	dummy_list = []
	fimg = pft.open(fits_name)
	prihdr = fimg[0].header
	scidata = fimg[0].data

	naxis2 = prihdr['NAXIS2']
	try:lampid = prihdr['LAMPID']
	except:lampid = 'NONE'
	try:grating = prihdr['GRATING']
	except:grating = 'NONE'
	try:gr_angle = prihdr['GR-ANGLE']
	except:gr_angle = 'NONE'
	try:object_n = prihdr['OBJECT']
	except:object_n = 'NONE'	
	dummy_list = [fits_name,naxis2,lampid,grating,gr_angle,object_n]
	return dummy_list

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

def trim(imname,x,y):
	fimg = pft.open(imname)
	prihdr = fimg[0].header
	scidata = fimg[0].data

	newscidata = scidata[y[0]-1:y[1],x[0]-1:x[1]]
	prihdr.update('NAXIS1', (x[1]-x[0]))
	prihdr.update('NAXIS2', (y[1]-y[0]))

	new_name = 'tr'+imname
	pft.writeto(new_name,data=newscidata,header=prihdr,clobber=True)
	os.system('mv %s history/' % (imname))

