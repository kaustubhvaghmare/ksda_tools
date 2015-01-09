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
#import pysynphot as ps
from RotationCurveLib import *


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

def SarithLimits(filelist):
	"""
	Input: List of Names of 1-d spectra.
	Output: w1, w2 the lower and upper wavelengths needed for calculation of sum.

	Program assumes equal length, does not check for the same.
	"""
	lefts = np.ones(len(filelist), float)
	rights = np.ones(len(filelist),float)
	for i,f in enumerate(filelist):
		hdulist = fits.open(f)
		wcs = pw.WCS(hdulist[0].header)
		lefts[i] = wcs.all_pix2world( [0], 1)[0][0] 
		rights[i] = wcs.all_pix2world( [len(hdulist[0].data)], 1)[0][0]

	edge1 = np.max(lefts)
	edge2 = np.min(rights)
	return edge1, edge2


############################
############################
def SumUpSpectra( spec_list, output_name, spec_err_list=False, output_err_filename=False, clean=True):
	"""
	Uses IRAF's sarith task to sum a list of spectra and generate an output spectrum
	Input: Python list of strings, each string being a file name.
		   Name of the output file to be generated.
		   clean=True, will delete old files are sum is generated.
	"""
	# We want to take care of the case where only one row is provided to us.
	if len(spec_list) == 1:
		os.rename(spec_list[0], output_name)
		if spec_err_list:
			os.rename(spec_err_list[0], output_err_filename)
		return

	# We want to estimate the wavelengths limits required.
	w1, w2 = SarithLimits(spec_list)
		
	# We want to start by adding two spectra first.
	iraf.noao.twodspec.longslit.sarith(input1=spec_list[0], input2=spec_list[1], op="+", output="sumspec_temp.fits", \
					clobber=True, format="onedspec")
	# Square the two error spectra
	if spec_err_list:
		#iraf.noao.twodspec.longslit.sarith(input1=spec_err_list[0], op="^",input2=2, output="temp_err1.fits", clobber=True, format="onedspec")
		#iraf.noao.twodspec.longslit.sarith(input1=spec_err_list[1], op="^",input2=2, output="temp_err2.fits", clobber=True, format="onedspec")
		iraf.stsdas.toolbox.imgtools.imcalc(input=spec_err_list[0], output="temp_err1.fits.0001.fits", equals="im1**2")
		iraf.stsdas.toolbox.imgtools.imcalc(input=spec_err_list[1], output="temp_err2.fits.0001.fits", equals="im1**2")
		#w1, w2 = SarithLimits("temp_err1.fits.0001.fits","temp_err2.fits.0001.fits") 
		iraf.noao.twodspec.longslit.sarith("temp_err1.fits.0001.fits", input2="temp_err2.fits.0001.fits", op="+", output="sumspec_temp_err.fits", \
					clobber=True, format="onedspec")
		os.system("rm temp_err1.fits.0001.fits temp_err2.fits.0001.fits")

	# What if there were no more than two spectra? The following loop will never execute.
	for spec in spec_list[2:]:
		#w1, w2 = SarithLimits("sumspec_temp.fits.0001.fits",spec)
		iraf.noao.twodspec.longslit.sarith("sumspec_temp.fits.0001.fits", input2=spec, op="+", output="sumspec_temp.fits", \
					clobber=True, format="onedspec")
		# Now, check for error frame corrections.
	if spec_err_list:
		for spec_err in spec_err_list:
			#iraf.noao.twodspec.longslit.sarith(input1=spec_err, op="^", input2=2, output="temp_err2.fits", clobber=True, format="onedspec")
			#iraf.images.imutil.imarith(operand1=spec_err, op="*", operand2=spec_err, result="temp_err2.fits.0001.fits")
			iraf.stsdas.toolbox.imgtools.imcalc(input=spec_err,output="temp_err2.fits.0001.fits", equals="im1**2")
			#w1,w2 = SarithLimits("sumspec_temp_err.fits.0001.fits","temp_err2.fits.0001.fits")
			iraf.noao.twodspec.longslit.sarith("sumspec_temp_err.fits.0001.fits", input2="temp_err2.fits.0001.fits", op="+",\
							output="sumspec_temp_err.fits", clobber=True, format="onedspec")
			os.system("rm temp_err2.fits.0001.fits")
		
	# Take error frame and divide by number of frames that went into it.
	iraf.noao.twodspec.longslit.sarith("sumspec_temp_err.fits.0001.fits", input2=len(spec_err_list), output="sumspec_temp_err.fits", op="/",\
						clobber=True, format="onedspec")
	# Take final error frame and square root it.
	#iraf.noao.twodspec.longslit.sarith(input1="sumspec_temp_err.fits.0001.fits", op="sqrt",output="sumspec_temp_err.fits", \
		#			clobber=True, format="onedspec")
	
	os.rename("sumspec_temp.fits.0001.fits", output_name)
	if spec_err_list:
		os.rename("sumspec_temp_err.fits.0001.fits", output_err_filename)
	# Perform clean-up unless user has turned that off!
	if clean:
		for s in spec_list:
			os.remove(s)
		for e in spec_err_list:
			os.remove(e)
	#return wl,wu


############################
############################
rotcurve = sys.argv[2]+"_fit.out"
rotcurve2 = sys.argv[2]+"_fit_z.out"

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
		for counter, row in enumerate(range(x1, x2)):
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
			redshift_at_aperture = splinez(row)
			# Deredshift the spectrum.
			iraf.noao.twodspec.longslit.dopcor(input=tempfile, output=tempfile, redshift = float(redshift_at_aperture))
			iraf.noao.twodspec.longslit.dopcor(input=temp_erfile, output=temp_erfile, redshift = float(redshift_at_aperture))
		# We now have several files present which have all been deredshifted.
		temp_aperture_files = [ "temp%d.fits" % i for i in range(len(range(x1,x2))) ]
		temp_aperture_errfiles = [ "temp%d_err.fits" % i for i in range(len(range(x1,x2))) ]
		aperture_filename = "%sWindow%d_Aperture%d.fits" % (filename_root, w, a) 
		aperture_err_filename = "%sWindow%d_Aperture%d_err.fits" % (filename_root, w, a) 
		SumUpSpectra( temp_aperture_files, aperture_filename, temp_aperture_errfiles, aperture_err_filename)
		print "Summing up all parts of aperture...done."
		position = (x1+x2)/2.0 - centroid
		aperture_map.write("%s\t%s\t%.2f\n" % (aperture_filename, aperture_err_filename, position))
