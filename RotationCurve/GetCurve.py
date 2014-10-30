import sys
import matplotlib.pyplot as plt
import numpy as np
import pyfits as pf
#from GetCurve_Config import *
from RotationCurveLib import *
import pywcs as pw
import scipy.optimize as so
from Fitting import *
from helcorr import *
from astropy.table import Table, Column
import os

# Parse command line arguments.
try:
	spectrum_file = sys.argv[1]
except:
	print("Insufficient number of arguments.")
	print("Usage: python GetCurve.py <2dspec_file>")
	sys.exit(2)

# Check if file exists, else throw error.
try:
	hdulist = pf.open(spectrum_file)
except:
	print "File does not exist or corrupt."
	sys.exit(3)

# Load the file, get data and extract wavelengths using WCS embedded info.
# Parse WCS information
header = hdulist[0].header
data = hdulist[0].data
wcs = pw.WCS(header)
wavelengths =  wcs.wcs_pix2sky( [ [i+0.5,1] for i in range( header["NAXIS1"] )] ,1 )
wavelengths = [ wavelengths[i][0] for i in range(len(wavelengths)) ]

# Determine the centroid row for displaying a spectrum.

centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )
central_row = int(round(centroid, 0))


# Display the file plot and enable input from the user.
fig1 = plt.figure(1)
plt.plot( wavelengths, data[central_row, :], color="blue")
plt.xlim( np.min(wavelengths), np.max(wavelengths) )
plt.title("Showing spectrum from row %d of %s" % (central_row, spectrum_file), fontsize=18)
plt.xlabel("Wavelengths [Angstroms]", fontsize=18)
plt.ylabel("Relative Flux", fontsize=18)
plt.grid()
c = Coordinate(fig1)
plt.show()
dumb_input = raw_input()
rest_wavelength = float(raw_input("Enter Rest Wavelength for Selected Feature: "))


x1, x2 = c.getcoords()
# To check if user entered coordinates in the right order.
if x1 > x2:
	x1,x2 = x2,x1

# Convert wavelength coordinates to pixels since we are dealing with those.
p1, p2 = Wavelength2Pixels(wcs, x1, x2)
p1, p2 = np.floor(p1), np.ceil(p2)


# Define some parameters. To be later moved into a separate config file.
stepsize = 5
nsum = 5
uppercut = 100
pixel_scale = 0.1 # arcsec per pixel	
speed_light = 3e5

# Create an array of pixels over which we will map the spectral feature.
# This, converted into wavelengths will serve as x in the Gaussian fit.
# Note that 0.5 addition in the xvar = .. step is because we want to take
# the wavelength for the mid-point of the pixel.
pixels = np.arange(p1, p2, 1)
xvar = Pixels2Wavelengths(wcs, *(pixels+0.5) )
xvar = np.array(xvar)

centroids = [] 
centroid_errs = []
spatial_points = []
for aperture in range(central_row, central_row + uppercut, stepsize):
	plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture )
	yvar = np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) 
	try:
		c, cerr = GaussFit(xvar, yvar)
		centroids.append(c)
		centroid_errs.append(cerr)
		spatial_points.append( aperture + float(stepsize)/2.0 )
	except:
		pass

for aperture in range(central_row, central_row - uppercut, -1*stepsize):
	plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture )
	yvar = np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) 
	try:
		c, cerr = GaussFit(xvar, yvar)
		centroids.append(c)
		centroid_errs.append(cerr)
		spatial_points.append( aperture + float(stepsize)/2.0 )
	except:
		pass

# Convert spatial row numbers to arc second distance from centre.
spatial_points_arc = (spatial_points - centroid)*pixel_scale

# Convert wavelength centroids into redshifts.
centroids_redshifts = (np.array(centroids) - rest_wavelength)/rest_wavelength
# Compute errors on redshifts.
centroids_redshifts_err = np.sqrt( (1.0/rest_wavelength)**2*np.array(centroid_errs)**2 )
# Convert wavelength centroids into km/s using rest-wavelength comparison. 
centroids_vel = (np.array(centroids) - rest_wavelength)/rest_wavelength * speed_light
# Compute heliocentric correction.
correction, hjd = helcorr(-20.810678, -32.376006, 1798, 4.269583, -55.780139, 2455084.537578)
# Apply heliocentric correction.
centroids_vel_helio = centroids_vel + correction
# Transform error bars accordingly.
centroids_vel_helio_err = np.sqrt( (speed_light/rest_wavelength)**2*np.array(centroid_errs)**2 + 1e-6)


# Plot rotation curve.
fig2 = plt.figure(2)
plt.errorbar( spatial_points_arc, centroids_vel_helio, centroids_vel_helio_err, fmt="o" )
plt.xlabel("Distance from Center [arcsec]", fontsize=18)
plt.ylabel("Heliocentric Radial Velocity [km/s]", fontsize=18)

# Output all rotation cruve data.
out = Table()
out.add_column( Column(spatial_points, name="Row") )
out.add_column( Column(spatial_points_arc, name="DFCG") )
out.add_column( Column(centroids, name="lambda") )
out.add_column( Column(centroid_errs, name="lambda_err") )
out.add_column( Column(centroids_vel, name="V") )
out.add_column( Column(centroids_vel_helio, name="Vhel") )
out.add_column( Column(centroids_vel_helio_err, name="Vhel_rr") )
out.add_column( Column(centroids_redshifts, name="z") )
out.add_column( Column(centroids_redshifts_err, name="z_err") )

# Ask user for file name to store rotation data.
filename = raw_input("Enter File Name to Save Data: ").rstrip()
out.write(filename, format="ascii", delimiter="\t")

# Display plot.
plt.grid()
fig2.show()
dum = raw_input()











