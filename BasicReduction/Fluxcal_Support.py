
# Supporter functions for the Fluxcal modules 
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

def apall(apname):
	outname = 'til'+apname
	refname = ""
	interap="yes"

	iraf.noao.twodspec.apextract.apall(input=apname,output=outname,apertur='',format='onedspec',referen=refname,profile='',
						interac=interap,find=interap,recente=interap,resize=interap,edit=interap,trace=interap,
						fittrac=interap,extract='yes',extras='no',review='no',line=500,nsum=30,
						lower=-5.,upper=12.,apidtab='',b_funct='chebyshev',b_order=1,b_sampl='-60:-40,40:60',
						b_naver=-3,b_niter=0,b_low_r=3.,b_high_=3.,b_grow=0.,width=5.,radius=10.,
						thresho=0.,nfind=1,minsep=35.,maxsep=100000.,order='increasing',aprecen='',
						npeaks='INDEF',shift='yes',llimit=-15,ulimit=15,ylevel=0.1,peak='yes',
						bkg='yes',r_grow=0.,avglimi='no',t_nsum=45,t_step=30,t_nlost=10,t_funct='spline3',
						t_order=3,t_sampl='*',t_naver=1,t_niter=0,t_low_r=3.,t_high_=3.,t_grow=0.,
						backgro='none',skybox=1,weights='none',pfit='fit1d',clean='no',saturat='INDEF',
						readnoi=0.,gain=1.,lsigma=4.,usigma=4.,nsubaps=1.,mode='ql')
	namesplit = string.split(outname,'.')
	os.system('ds9 %s &' % (namesplit[0]+'.0001.fits'))
	satis = input_str("Are satisfied with the extracted spectra (y|n)?")
	os.system("kill -9 `pidof ds9`")
	if (satis == 'n'):
		askdel = input_str("Delete the spectra (y|n)?")
		if (askdel == 'y'):
			os.system('rm %s' % (namesplit[0]+'.0001.fits'))
	else:
		os.system('mv %s history/' % (apname))
	return


## The standardisation module.
def do_calibration(std_name,starname):
	outname = 'til'+std_name
	namesplit = string.split(outname,'.')
	std_specname = namesplit[0]+'.0001.fits'
	iraf.noao.twodspec.longslit.standard(input=std_specname, output="std", samestar="yes", beam_switch="no", apertures="", bandwidth="INDEF",
					    bandsep="INDEF", fnuzero=3.68e-20, extinction=irafhome+"suth_extinct.dat",interact="yes", 
					    star_name=starname)
	
	iraf.noao.twodspec.longslit.sensfunc(standards="std", sensitivity="sens", apertures="", ignoreaps="yes", logfile="logfile", 
					        extinction=irafhome+"suth_extinct.dat", function="spline3", order=12, interactive="yes")

						
# The calibration module.
def calibrate(science):
	outname = "cal_" + science

	iraf.noao.twodspec.longslit.calibrate(input=science, output=outname, extinct="yes", flux="yes", extinction=irafhome+"suth_extinct.dat", 
						ignoreaps="yes", sensitivity="sens", fnu="no")
	
	
					    








