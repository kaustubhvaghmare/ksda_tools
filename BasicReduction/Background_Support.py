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

def background(sciname,filename):
	print "Working on %s." % sciname
	backn = 'b'+sciname
	os.system("ds9 %s -zscale &" % (sciname))
	iraf.noao.twodspec.longslit.background(input=sciname,output=backn,axis=2,interac='yes',sample='*',naverag=1,
					functio='spline3',order=3,low_rej=2.,high_re=1.5,niterat=1,grow=0.,
					graphic='stdgraph',cursor='',mode='al')
	try:
		os.system("kill -9 `pidof ds9`")
	except:
		pass
	print "Please review background subtracted image."
	os.system('ds9 %s -zscale &' % (backn))
	os.system('mv %s history/' % (sciname))
	raw_input("Press any key to proceed to next image.")
	try:
		os.system("kill -9 `pidof ds9`")
	except:
		pass
	return

