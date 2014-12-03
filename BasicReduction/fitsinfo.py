#! /usr/bin/python
''' Mini code to output the properties of fits files inside the folder you are working '''


import numpy as np
import pyfits as pft
import os,sys,glob,string
import math as mth

#function to read the header of fits files and output a few critical info needed in the program
def info_fits(fits_name):
	dummy_list = []
	fimg = pft.open(fits_name)
	prihdr = fimg[0].header
	try:prihdr1 = fimg[1].header
	except:prihdr1 = 'none'
	try:naxis2 = prihdr1['NAXIS2']
	except: naxis2 = 'none'
	try:naxis2 = prihdr['NAXIS2']
	except: naxis2 = 'none'
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

files = glob.glob('*fits')
#renaming the files to more appropriate names!
print ' name \t naxis2 \t lampid \t grating \t gr_angle \t object_n'
for i in range(0,len(files)):
	print info_fits(files[i])
