"""
A collection of functions to be used by GetCurve2.py in order to 
perform various operations on it.
"""

from astropy.io import fits
import astropy.wcs as pw
import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import GaussianModel, ConstantModel
import string
import ConfigParser
from astropy import time

# Define a simple function to read a 2-d spectrum.
def read_spectrum(filename):
	"""
	Input: Name of the file.
	Output: header object, data object, wavelengths for pixel centres.
		In case the loading of the file fails, you get a -1,-1, -1 tuple.
	"""
	try:
		hdulist = fits.open(filename)
	except:
		return -1,-1,-1
	header = hdulist[0].header
	data = hdulist[0].data
	wcs = pw.WCS(header)
	wavelengths =  wcs.wcs_pix2world( [ [i+0.5,1] for i in range( header["NAXIS1"] )] ,1 )
	wavelengths = [ wavelengths[i][0] for i in range(len(wavelengths)) ]

	return header, data, wavelengths

# A Matplotlib Connection Event Object to query the coordinates of plots.
class Coordinate:
	def __init__(self, figure):
		self.x = []
		self.figure = figure
		figure.canvas.mpl_connect("key_press_event", self)

	def __call__(self, event):
		if len(self.x) < 2 and event.key == "m":
			self.x.append( event.xdata )
			print "Coordinate %d added = %.2f" % (len(self.x), event.xdata)
		if len(self.x) == 2:
			print "These coordinates will be used for Rotation Curve. Press any key to continue."
			plt.close(self.figure)

	def getcoords(self):
		if self.x[0] < self.x[1]:
			return self.x
		else:
			return self.x[::-1]


# Define a centroiding algorithm, use it get a central row.
def GetCentroid( spatial_points, intensities):
	return np.sum(spatial_points * intensities)/np.sum(intensities)

def Pixels2Wavelengths(header, *args):
	wcs = pw.WCS(header)
	pixels = [ [args[i]+0.5,1] for i in range(len(args)) ]
	coords = wcs.wcs_pix2world( pixels, 1)
	coords = [ coords[i][0] for i in range(len(pixels)) ]
	return coords

def Wavelength2Pixels(header, *args):
	wcs = pw.WCS(header)
	coords = [ [args[i],1] for i in range(len(args)) ]
	pixels = wcs.wcs_world2pix( coords, 1)
	pixels = [ pixels[i][0] for i in range(len(pixels)) ]
	return pixels

def FitGaussians(x, y, y_err, num):
	"""
	Input: x, generally wavelengths.
	       y, generally flux
		   y_err, errors on y
		   num: number of Gaussians to fit
	Output: Mean Centroid, Error on Mean Centroid
	"""
	# First, define a model.
	mod = GaussianModel(prefix="g1_")
	components = 1
	while components != num:
		mod = mod + GaussianModel(prefix="g%d_" % (components+1))
		components += 1
	mod = mod + ConstantModel()
	pars = mod.make_params()

	# Some initial estimates.
	bins_on_x = np.linspace(x.min(),x.max(),num+1)
	x_digitized = np.digitize(x,bins_on_x)
	
	for i in range(num):
		prefix = "g%d_" % (i+1)
		guess_center = np.mean(x[x_digitized == (i+1)])
		pars[prefix+"center"].set(guess_center, min=x.min(), max=x.max())
		pars[prefix+"amplitude"].set(np.sum(y)/num)
		pars[prefix+"sigma"].set(5)
		pars["c"].set(np.mean(y))

	# Great, model has been defined, parameters defined. 
	results = mod.fit(y, pars, x=x, weights=1/y_err)
	params = results.params
	centers = np.ones(num, dtype=float)
	centers_err = np.ones(num, dtype=float)
	for i in range(num):
		param_name = "g%d_center" % (i+1)
		centers[i] = params[param_name].value
		centers_err[i] = params[param_name].stderr
	
	ind_sort = np.argsort(centers)
	return centers, centers_err

# A basic y/n questioner.
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

# An object input file parser.
def ObjectInformation(header):
	h,m,s = [ float(i) for i in header["RA"].split(":") ]
	ra = 15*(h + m/60.0 + s/3600.0)
	d,m,s = [ float(i) for i in header["DEC"].split(":") ]
	if d<0:
		dec = (d + m/60.0 + s/3600.0)*-1.0
	else:
		dec = (d + m/60.0 + s/3600.0)
	
	time_string = header["DATE-OBS"] + " " + header["UTC-OBS"]
	t = time.Time(time_string)
	jd = t.jd

	return ra, dec, jd
	


# Observatory Information
def ObservatoryInformation(filename):
	c = ConfigParser.ConfigParser()
	c.read(filename)
	longi = c.getfloat("ObservatoryInformation", "longitude")
	lat = c.getfloat("ObservatoryInformation", "latitude")
	alt = c.getfloat("ObservatoryInformation", "altitude")
	
	return longi, lat, alt

# Curve determination parameters.
def CurveParams(filename):
	c = ConfigParser.ConfigParser()
	c.read(filename)
	
	stepsize = c.getint("CurveParameters","stepsize")
	nsum = c.getint("CurveParameters","nsum")
	uppercut = c.getint("CurveParameters", "uppercut")
	pixel_scale = c.getfloat("CurveParameters", "pixelscale")
	speed_light = c.getfloat("CurveParameters","speedlight")

	return stepsize, nsum, uppercut, pixel_scale, speed_light
