"""
A collection of functions to be used by GetCurve2.py in order to 
perform various operations on it.
"""

from astropy.io import fits
import astropy.wcs as pw
import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import GaussianModel, LinearModel
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
	wavelengths = np.array( [ wavelengths[i][0] for i in range(len(wavelengths)) ] )

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

def FitGaussians(x, y, y_err, num, rest_wavelengths=False, mode="centroid"):
	"""
	Input: x, generally wavelengths.
	       y, generally flux
		   y_err, errors on y
		   num: number of Gaussians to fit
		   rest_wavelengths: used for computing constraints on Gaussian fits
	Output: Mean Centroid, Error on Mean Centroid if used in default,
			if mode = "width", will return computed sigma and sigma_err.
	"""
	# First, define a model.
	mod = GaussianModel(prefix="g1_")
	components = 1
	while components != num:
		mod = mod + GaussianModel(prefix="g%d_" % (components+1))
		components += 1
	mod = mod + LinearModel()
	pars = mod.make_params()

	# Some initial estimates.
	bins_on_x = np.linspace(x.min(),x.max(),num+1)
	x_digitized = np.digitize(x,bins_on_x)
	
	for i in range(num):
		prefix = "g%d_" % (i+1)
		guess_center = np.mean(x[x_digitized == (i+1)])
		if i == 0:
			pars[prefix+"center"].set(guess_center, min=x.min(), max=x.max())
			pars[prefix+"amplitude"].set(np.sum(y)/num)
			pars[prefix+"sigma"].set(5,min=0)
		else:
 			mu_del = rest_wavelengths[i]- rest_wavelengths[0] # difference
			pars[prefix+"center"].set(guess_center, min=x.min(), max=x.max(), expr='g1_center + %f' % mu_del )
			pars[prefix+"amplitude"].set(np.sum(y)/num)
			pars[prefix+"sigma"].set(5, expr='g1_sigma')

	pars["intercept"].set(np.mean(y))
	pars['slope'].set(0.1)

	# Great, model has been defined, parameters defined. 
	results = mod.fit(y, pars, x=x, weights=1/y_err)
	params = results.params
	centers = np.ones(num, dtype=float)
	centers_err = np.ones(num, dtype=float)
	widths = np.ones(num, dtype=float)
	widths_err = np.ones(num, dtype=float)
	for i in range(num):
		param_name1 = "g%d_center" % (i+1)
		param_name2 = "g%d_sigma" % (i+1)
		centers[i] = params[param_name1].value
		centers_err[i] = params[param_name1].stderr
		widths[i] = params[param_name2].value
		widths_err[i] = params[param_name2].stderr

	ind_sort = np.argsort(centers)

	if mode=="centroid":
		return centers, centers_err
	elif mode=="width":
		return widths, widths_err
	elif mode=="both":
		return centers, centers_err, widths, widths_err
	else:
		print("Warning: Invalid mode argument. Returning centroids.")
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
	try:
		h,m,s = [ float(i) for i in header["RA"].split(":") ]
		ra = 15*(h + m/60.0 + s/3600.0)
	except:
		ra = float(header["RA"])
	try:
		d,m,s = [ float(i) for i in header["DEC"].split(":") ]
		if d<0:
			dec = (d + m/60.0 + s/3600.0)*-1.0
		else:
			dec = (d + m/60.0 + s/3600.0)
	except:
		dec = float(header["DEC"])
	
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

def FilterErrorBars(x_err):
	"""
	Input: x and x_err
	Output: x and x_err filtered of insane error bars.
	Notes: Not meant to be a general function. Have some specific nuances
	meant for this.
	"""
	x_err = np.array(x_err)
	med = np.median(x_err)
	desired = ( x_err <= 2*med )
	return desired


# Function for extraction of 1-d spectrum from a 2-d spectrum.
def Extract1d(spec_file, row_num, col1, col2, out_file, clobber=True):
	"""
	Input: The name of the input 2d spectrum, the row of extraction and the output file.
	Output: None. Just creates the relevant file.
	"""
	if col1 > col2:
		print("Col1 cannot be greater than Col2")
		return 
	hdulist = fits.open(spec_file)
	data = hdulist[0].data[row_num,col1:col2+1 ]
	header = hdulist[0].header
	for key in ["CTYPE2", "CDELT2", "CD2_2", "LTM2_2", "WAT2_001", "DISPAXIS"]:
		header.remove(key)
	header["WCSDIM"] = 1
	header["WAT0_001"] = "system=equispec"
	header["CRPIX1"] = col1 - 1 - header["CRPIX1"]
	header["LTV1"] = 1 - col1

	fits.writeto(out_file, data, header, clobber=clobber)

def DopplerCorrect(spec_file, redshift, out_file, clobber=True):
	"""
	Input: The input 1d spec_file to be Doppler Corrected.
		   The redshift to be applied.
		   The output file name.
	Output: None. Just creates a new file with Doppler correction applied.
	"""
	hdulist = fits.open(spec_file)
	data = hdulist[0].data
	header = hdulist[0].header
	header["CRVAL1"] = header["CRVAL1"] / (1.0+redshift)
	header["CDELT1"] = header["CDELT1"] / (1.0+redshift)
	header["CD1_1"] = header["CD1_1"] / (1.0+redshift)

	fits.writeto(out_file, data, header, clobber=clobber)

def MeanSigmaClipper(x, sigma=3, compliment=False,niter=10):
	"""
	Input: x, the array on which simple sigma clipping needs to be done.
		   sigma, the level of clipping
		   compliment, by default those that remain post rejection, their
		   indices are returned, but this was cause indices of rejected
		   to be returned.
	Output: Indices of those elements in x that withstand rejection or
		   otherwise if the compliment option is selected.
	"""
	accepted = np.ones( len(x), dtype=bool)

	for i in range(niter):
		xt = x[accepted]
		mean = np.mean(xt)
		std = np.std(xt)
		filt = (x <= (mean + sigma*std) ) & (x >= (mean - sigma*std) )
		accepted = np.bitwise_and( accepted, filt )
	
	if not compliment:
		return accepted
	else:
		return np.bitwise_not(accepted)
	
