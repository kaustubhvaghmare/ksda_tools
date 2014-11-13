"""
This program will do the following.
Mark apertures on a given 2d spectrum
For each aperture, 
	* use the obtained solution for its rotation curve to derotate
	* sum up derotated aperture rows.
	* dump them out as a 1d fits file.
"""

import sys
import astropy.io.fits as fits
import astropy.wcs as wcs
import pickle
import os
import ApertureLibrary as al
from pyraf import iraf
import pysynphot as ps
import numpy as np

# Start by validating inputs.
if len(sys.argv) != 3:
	print("Incorrect Usage.")
	print("Usage is: python DerotExtract.py <2dspecfile> <rotationcurve>")
	sys.exit(2)

# Check if files specified are genuine.
filename = sys.argv[1]
if not os.path.isfile(filename):
	print("2d spectrum FITS file does not exist.")
	sys.exit(3)
filename_root = filename.replace(".fits", "")

############################
############################
rotcurve = sys.argv[2]+"_spline.out"
rotcurve2 = sys.argv[2]+"_spline_z.out"

if not os.path.isfile(rotcurve) or not os.path.isfile(rotcurve2):
	print("Fitting has not been carried for this curve. Please fit the solution first using FitRotationCurve.py")
	sys.exit(4)

# Load the solution.
splinez = pickle.load(open(rotcurve2))

# Take the user through the following steps
# Display spatial profile.
# Mark windows / apertures
# Get aperture coordinate array.
aperture_coordinates = al.get_aperture(filename)

# aperture_coordinates[0] = array of apertures for first window.
# This is how aperture_coordinates is structured.
(no_windows, no_apertures) = aperture_coordinates.shape

hdulist = fits.open(filename)
data = hdulist[0].data
header = hdulist[0].header

wcs_object = wcs.WCS(header)
pixels = np.array( [ [i+0.5,1] for i in range(data.shape[1]) ] )
wavelengths = wcs_object.all_pix2sky(pixels, 0).T[0]


aperture_map = open(filename_root+"_aper_map.dat", "w")

for w in range(no_windows):
	print "Working on Window %d:" % w
	for a in range(no_apertures-1):
		x1, x2 = int(aperture_coordinates[w, a]), int(aperture_coordinates[w, a+1])
		print "\tWorking on Aperture %d defined from %d to %d." %(a,x1,x2)
		summed_spectrum = ps.ArraySpectrum(flux=np.zeros(data.shape[1]), wave=wavelengths)

		for counter, row in enumerate(range(x1, x2)):
			# Get data slice.
			spec_data = data[row,:]
			# Make a spectrum out of it using original header.
			spectrum = 	ps.ArraySpectrum(flux=spec_data, wave=wavelengths)
			# Use Spline to get wavelength associated with 
			redshift_at_aperture = splinez(row)
			# Deredshift the spectrum.
			spectrum.redshift(-1*redshift_at_aperture)
			summed_spectrum = spectrum + summed_spectrum
			
		aperture_filename = "%sWindow%d_Aperture%d.fits" % (filename_root, w, a) 
		summed_spectrum.writefits(aperture_filename)
		print "Summing up all parts of aperture...done."
		aperture_map.write("%s\t%.2f\n" % (aperture_filename, (x1+x2)/2.0) )


	


		

	



