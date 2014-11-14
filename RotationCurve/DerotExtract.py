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
def SumUpSpectra( spec_list, output_name, clean=True):
	"""
	Uses IRAF's sarith task to sum a list of spectra and generate an output spectrum
	Input: Python list of strings, each string being a file name.
		   Name of the output file to be generated.
		   clean=True, will delete old files are sum is generated.
	"""
	# We want to take care of the case where only one row is provided to us.
	if len(spec_list) == 1:
		os.rename(spec_list[0], output_name)
		return
	# We want to start by adding two spectra first.
	iraf.noao.twodspec.longslit.sarith(input1=spec_list[0], input2=spec_list[1], op="+", output="sumspec_temp.fits", \
					clobber=True, format="onedspec")

	# What is there were no more than two spectra? The following loop will never execute.
	for spec in spec_list[2:]:
		iraf.noao.twodspec.longslit.sarith("sumspec_temp.fits.0001.fits", input2=spec, op="+", output="sumspec_temp.fits", \
					clobber=True, format="onedspec")

	os.rename("sumspec_temp.fits.0001.fits", output_name)
	# Perform clean-up unless user has turned that off!
	if clean:
		for s in spec_list:
			os.remove(s)


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
no_windows = aperture_coordinates.shape[0]
no_apertures = (aperture_coordinates[0] != 0).sum()

hdulist = pf.open(filename)
data = hdulist[0].data
header = hdulist[0].header

# Determine centroid of the spatial profile for the 2d spectra.
centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )


aperture_map = open(filename_root+"_aper_map.dat", "w")

for w in range(no_windows):
	print "Working on Window %d:" % w
	for a in range(no_apertures-1):
		x1, x2 = int(aperture_coordinates[w, a]), int(aperture_coordinates[w, a+1])
		print "\tWorking on Aperture %d defined from %d to %d." %(a+1,x1,x2)
		for counter, row in enumerate(range(x1, x2)):
			# Get data slice.
			spec_data = data[row,:]
			# Make a spectrum out of it using original header.
			tempfile = "temp%d.fits" % counter
			pf.writeto(tempfile, data = spec_data, header=header)
			# Use Spline to get wavelength associated with 
			redshift_at_aperture = splinez(row)
			# Deredshift the spectrum.
			iraf.noao.twodspec.longslit.dopcor(input=tempfile, output=tempfile, redshift = float(redshift_at_aperture))
		# We now have several files present which have all been deredshifted.
		temp_aperture_files = [ "temp%d.fits" % i for i in range(len(range(x1,x2))) ]
		aperture_filename = "%sWindow%d_Aperture%d.fits" % (filename_root, w, a) 
		SumUpSpectra( temp_aperture_files, aperture_filename)
		print "Summing up all parts of aperture...done."
		position = (x1+x2)/2.0 - centroid
		aperture_map.write("%s\t%.2f\n" % (aperture_filename, position))
