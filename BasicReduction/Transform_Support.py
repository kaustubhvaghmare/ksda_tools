# Support files to the Transform Module.

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

def transform(sciname,filename):
	# taking the arc file name and using it as the reference input for transformation
	namesplit = string.split(filename,'.')
	trans = 't'+sciname
	iraf.noao.twodspec.longslit.transform(input=sciname,output=trans,minput='',moutput='',fitnames=namesplit[0],
					databas='database',
					interpt='linear',x1='INDEF',x2='INDEF',dx='INDEF',nx='INDEF',xlog='no',y1='INDEF',
					y2='INDEF',dy='INDEF',ny='INDEF',ylog='no',flux='yes',blank='INDEF',
					logfile='STDOUT,logfile',mode='ql')
	os.system('ds9 %s -zscale &' % (trans))
	os.system('mv %s history/' % (sciname))



