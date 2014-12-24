# Import libraries.
import numpy as np
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import stsdas

import pyfits as pft
import os,sys,glob,string,time
import math as mth
import pickle

# Define the IRAF home directory.
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

# Define the combine function.
def Combine(list_string, com_name):
	no_images = len(list_string.split(","))
	iraf.images.immatch.imcombine(input=list_string, output=com_name, sigma=com_name+"_tempstd.fits", combine="average", 
							reject="avsigclip", project="no", outtype="real", lsigma=3.0, hsigma=3.0)
	# Modify _tempstd.fits image to truly reflect the error bars.
	iraf.stsdas.toolbox.imgtools.imcalc(input=com_name+"_tempstd.fits", output=com_name+"_tempstd2.fits", equals="im1**2")
	iraf.stsdas.toolbox.imgtools.imcalc(input=com_name+"_tempstd2.fits", output=com_name+"_tempstd3.fits", equals="im1/%d" % no_images)
	iraf.stsdas.toolbox.imgtools.imcalc(input=com_name+"_tempstd3.fits", output=com_name+"_err.fits", equals="sqrt(im1)")
	os.system("rm %s_tempstd*.fits" % com_name)
