"""
This program will do the following.
Mark apertures on a given 2d spectrum
For each aperture, 
	* use the obtained solution for its rotation curve to derotate
	* sum up derotated aperture rows.
	* dump them out as a 1d fits file.
"""

import sys
import pyfits as pf
import pickle
import os
import ApertureLibrary as al
from pyraf import iraf
import pysynphot as ps
from RotationCurveLib import *
from DerotLibrary import *
from numpy.polynomial.chebyshev import chebval

# Start by validating inputs.
if len(sys.argv) != 4:
	print("Incorrect Usage.")
	print("Usage is: python DerotExtract.py <2dspecfile> <rotationcurve> <2dspec_errorfile>")
	sys.exit(2)

# Check if files specified are genuine.
filename = sys.argv[1]
if not os.path.isfile(filename):
	print("2d spectrum FITS file does not exist.")
	sys.exit(3)
filename_root = filename.replace(".fits", "")

err_filename = sys.argv[3]
if not os.path.isfile(filename):
	print("2d spectrum FITS file does not exist.")
	sys.exit(3)
errfilename_root = filename.replace(".fits", "")

############################
############################

rotcurve = sys.argv[2]+"_fit.out"
rotcurve2 = sys.argv[2]+"_fit_z.out"

if not os.path.isfile(rotcurve) or not os.path.isfile(rotcurve2):
	print("Fitting has not been carried for this curve. Please fit the solution first using FitRotationCurve.py")
	sys.exit(4)

# Load the solution.
interp_coeff_z = pickle.load(open(rotcurve2))

# Take the user through the following steps
# Display spatial profile.
# Mark windows / apertures
# Get aperture coordinate array.
aperture_coordinates = al.get_aperture(filename)

# aperture_coordinates[0] = array of apertures for first window.
# This is how aperture_coordinates is structured.
no_windows = aperture_coordinates.shape[0]
no_apertures = (aperture_coordinates[0] != 0).sum()

hdulist = pf.open(filename)
data = hdulist[0].data
header = hdulist[0].header

er_hdulist = pf.open(err_filename)
er_data = er_hdulist[0].data
er_header = er_hdulist[0].header

# Determine centroid of the spatial profile for the 2d spectra.
centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )


aperture_map = open(filename_root+"_aper_map.dat", "w")

for w in range(no_windows):
	print "Working on Window %d:" % w
	for a in range(no_apertures-1):
		x1, x2 = int(aperture_coordinates[w, a]), int(aperture_coordinates[w, a+1])
		print "\tWorking on Aperture %d defined from %d to %d." %(a+1,x1,x2)
		for counter, row in enumerate(range(x1, x2+1)):
			# Get data slice.
			spec_data = data[row,:]
			spec_err = er_data[row,:]
			# Make a spectrum out of it using original header.
			tempfile = "temp%d.fits" % counter
			temp_erfile = "temp%d_err.fits" % counter
			
			iraf.noao.twodspec.longslit.scopy(input = filename+"[*,%d]" % row, output=tempfile, format="onedspec")
			os.system("mv %s %s" % (tempfile+".0001.fits", tempfile))
			iraf.noao.twodspec.longslit.scopy(input = err_filename+"[*,%d]" % row, output=temp_erfile, format="onedspec")
			os.system("mv %s %s" % (temp_erfile+".0001.fits", temp_erfile))
#			pf.writeto(tempfile, data = spec_data, header=header)
#			pf.writeto(temp_erfile, data = spec_err, header=er_header)
			# Use Spline to get wavelength associated with 
			redshift_at_aperture = chebval(row, interp_coeff_z)
			# Deredshift the spectrum.
			iraf.noao.twodspec.longslit.dopcor(input=tempfile, output=tempfile, redshift = float(redshift_at_aperture))
			iraf.noao.twodspec.longslit.dopcor(input=temp_erfile, output=temp_erfile, redshift = float(redshift_at_aperture))
		# We now have several files present which have all been deredshifted.
		temp_aperture_files = [ "temp%d.fits" % i for i in range(len(range(x1,x2+1))) ]
		temp_aperture_errfiles = [ "temp%d_err.fits" % i for i in range(len(range(x1,x2+1))) ]
		aperture_filename = "%sWindow%d_Aperture%d.dat" % (filename_root, w, a) 
		aperture_err_filename = "%sWindow%d_Aperture%d_err.dat" % (filename_root, w, a) 
		SumUpSpectra( temp_aperture_files, aperture_filename, temp_aperture_errfiles, aperture_err_filename)
		print "Summing up all parts of aperture...done."
		position = (x1+x2)/2.0 - centroid
		aperture_map.write("%s\t%s\t%.2f\n" % (aperture_filename, aperture_err_filename, position))

os.system("rm temp*.fits")
