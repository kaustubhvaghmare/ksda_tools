# Support code fof Identify.py

# Basic library imports.

import numpy as np
from pyraf import iraf
from pyraf.iraf import noao
from pyraf.iraf import stsdas

import pyfits as pft
import os,sys,glob,string,time
import math as mth
import pickle

from IRAFHome import *

def identify(properties,index,irafhome):
	satis = 'n'
	while satis == 'n':
	        properties = np.array(properties)
		print properties
		filename = properties[index,0]
		lampid = properties[index,2]
		path = irafhome+lampid+'.txt'
	 	print path,filename
		iraf.noao.twodspec.longslit.identify(images=filename, section='middle line', databas='database',
			            	coordli=path, units='',nsum=10,match=-3., maxfeat=50, zwidth=100.,ftype='emission',
				    	fwidth=5., cradius=6., thresho=0., minsep=2.,functio='spline3',order=3,sample= '*',
				    	niterat=0, low_rej=3., high_re=3., grow=0.,autowri='no', graphic='stdgraph', cursor='',  	       					    	crval='',cdelt='',aidpars='',mode='ql')
		iraf.noao.twodspec.longslit.reidentify(referenc=filename,images=filename,interac='no',section='middle line',newaps='yes',
					overrid='no',refit='no',trace='no',step=10.,nsum=10.,shift=0.,search=5.,nlost=10.,
					cradius=5.,thresho=0.,addfeat='no',coordli=path,match=-3,maxfeat=50,minsep=2,
					databas='database',logfile='logfile',plotfil='',verbose='yes',graphic='stdgraph',
					cursor='',answer='yes',crval='',cdelt='',aidpars='',mode='ql')
		namesplit = string.split(filename,'.')
		iraf.noao.twodspec.longslit.fitcoords(images=str(namesplit[0]),fitname='',interac='yes',combine='no',
					databas='database',deletio='deletions.db',
					functio='chebyshev',xorder=6,yorder=6,logfile='STDOUT,logfile',plotfil='plotfile',
					graphic='stdgraph',cursor ='',mode='ql')
		trans = 't'+filename
		iraf.noao.twodspec.longslit.transform(input=filename,output=trans,minput='',moutput='',fitnames=namesplit[0],
					databas='database',
					interpt='spline3',x1='INDEF',x2='INDEF',dx='INDEF',nx='INDEF',xlog='no',y1='INDEF',
					y2='INDEF',dy='INDEF',ny='INDEF',ylog='no',flux='yes',blank='INDEF',
					logfile='STDOUT,logfile',mode='ql')
		os.system('ds9 %s -zscale &' % (trans))
		satis = input_str("Are satisfied with the transformed spectra (y|n)? :")
		if (satis == 'n'):
			askdel = input_str("Delete the spectra and database (y|n)? :")
			if (askdel == 'y'):
				os.system('rm -r database/ %s deletion.db' % (trans))
		os.system('mv %s %s history/' % (filename,trans))
		try:
			os.system('kill -9 `pidof ds9`')
		except:
			pass
	return

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

