import numpy as np
import matplotlib.pyplot as plt
from Fitting import *


# Define a centroiding algorithm, use it get a central row.
def GetCentroid( spatial_points, intensities):
	return np.sum(spatial_points * intensities)/np.sum(intensities)

# # Define an event for accepting input coordinates via clicking.
class Coordinate:
	def __init__(self, figure):
		self.x = []
		self.figure = figure
		figure.canvas.mpl_connect("button_press_event", self)

	def __call__(self, event):
		if len(self.x) < 2:
			self.x.append( event.xdata )
			print "Coordinate %d added = %.2f" % (len(self.x), event.xdata)
		if len(self.x) == 2:
			print "These coordinates will be used for Rotation Curve. Press any key to continue."
			plt.close(self.figure)

	def getcoords(self):
		return self.x

# Define simple functions to convert pixels to wavelength and vice versa.
def Wavelength2Pixels(wcs, *args):
	coords = [ [args[i],1] for i in range(len(args)) ]
	pixels = wcs.wcs_sky2pix( coords, 1)
	pixels = [ pixels[i][0] for i in range(len(pixels)) ]
	return pixels

def Pixels2Wavelengths(wcs, *args):
	pixels = [ [args[i]+0.5,1] for i in range(len(args)) ]
	coords = wcs.wcs_pix2sky( pixels, 1)
	coords = [ coords[i][0] for i in range(len(pixels)) ]
	return coords

# We will fit profiles with Gaussian, so define a fitting function for that.
def GaussFit(x, y, yerr=None, p0=None):
	def gaussian(x, a, mu, sigma, offset):
		return a*np.exp((x-mu)**2/(-2*sigma**2) ) + offset

	if p0 is None:
		p0 = [ np.sum(y), np.mean(x), np.std(x), np.max(y)-np.min(y) ]

	gau_params, gau_unc, rchi2, dof= general_fit( gaussian, x, y, p0, yerr)
	return gau_params[1], gau_unc[1]

