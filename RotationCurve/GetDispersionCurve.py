"""
Program to determine velocity dispersion curve.
Usage: python GetDispersionCurve.py <2dspec_file> <rot_curve> <2dspecerr_file>
"""

import sys
import matplotlib.pyplot as plt
from MatplotlibCustom import *
from RotationCurveLib import *
from astropy.table import Table, Column
import pickle
from pyraf import iraf
import os

######################################
# Parse command line arguments.
######################################
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

rotcurve = sys.argv[2]+"_fit.out"
rotcurve2 = sys.argv[2]+"_fit_z.out"

if not os.path.isfile(rotcurve) or not os.path.isfile(rotcurve2):
	print("Fitting has not been carried for this curve. Please fit the solution first using FitRotationCurve.py")
	sys.exit(4)

#######################################################
# Load Files and report any errors found in the process.
#######################################################
header, data, wavelengths = read_spectrum( filename )
if header == -1:
	print("File corrupt or does not exist.")
	sys.exit(3)
er_header, er_data, er_wavelengths = read_spectrum( err_filename )
if er_header == -1:
	print("Error file corrupt or does not exist.")
	sys.exit(4)
splinez = pickle.load(open(rotcurve2))

########################################################################
# Displaying the spectrum and allowing user to choose features to model.
########################################################################

dopcor = input_str("Do you want to display spectrum in an approximately deredshifted frame? ")
if dopcor == "y":
	redshift = float(raw_input("Approximate redshift, please? "))
else:
	redshift = 0.0

# Determine the centroid row for displaying a spectrum.
centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )
central_row = int(round(centroid, 0))

# Enclosing entire process of using spectral lines in a while loop
# which will be executed indefinity till the user says "n" to the 
# question - "Do you want another rotation curve?"

while True:
	# Display the file plot and enable input from the user.
	fig1 = plt.figure(1)
	plt.plot( wavelengths/(1+redshift), data[central_row, :], color="blue")
	plt.xlim( np.min(wavelengths), np.max(wavelengths) )
	plt.title("Showing spectrum from row %d of %s" % (central_row, filename), fontsize=18)
	plt.xlabel("Wavelengths [Angstroms]", fontsize=18)
	plt.ylabel("Relative Flux", fontsize=18)
	plt.grid()
	c = Coordinate(fig1)
	print("Went past c = Coordinate(fig1)");
	plt.show()
	print("Went pas show()")
	dumb_input = raw_input()
	print("Got dumb input")
	x1, x2 = c.getcoords()
	x1 = x1*(1+redshift)
	x2 = x2*(1+redshift)
	
	# Prompt the user for number of Gaussian components needed to model the selected spectrum.
	no_components = int(raw_input("Enter number of Gaussian Components to model: "))
	
	# Prompt the user for rest wavelengths.
	rest_wavelengths = np.ones(no_components)
	for i in range(no_components):
		rest_wavelengths[i] = float(raw_input("Enter Rest Wavelength for Selected Feature(s): "))

	########################################################
	# Rotation Curve Extraction
	########################################################
	
	# First off, we need to access data from the spectrum but for user
	# convenience we have dealt with wavelengths. Let's get the pixels.
	p1, p2 = Wavelength2Pixels(header, x1, x2)
	p1, p2 = np.floor(p1), np.ceil(p2)
	
	# Read configuration file.
	stepsize, nsum, uppercut, pixel_scale, speed_light = CurveParams(config_file)
	
	# Create an array of pixels over which we will map the spectral feature.
	# This, converted into wavelengths will serve as x in the Gaussian fit.
	# Note that 0.5 addition in the xvar = .. step is because we want to take
	# the wavelength for the mid-point of the pixel.
	pixels = np.arange(p1, p2, 1)
	xvar = Pixels2Wavelengths(header, *(pixels+0.5) )
	xvar = np.array(xvar)  # yvar will be determined for each aperture.

	widths = []
	widths_err = []
	centroids = []
	centroids_err = []
	spatial_points = []

	# Move over each aperture.
	for aperture in range(central_row, central_row + uppercut, stepsize):
		plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture ) # Diagnostic
		
		# Not to construct yvar for this aperture. This involves three steps.
		for row in range(aperture, aperture+nsum+1):
			
			# 1. Take each row and throw in a separate fits file.
			tempfile = "temp%d.fits" % counter
			temp_erfile = "temp%d_err.fits" % counter

			iraf.noao.twodspec.longslit.scopy(input = filename+"[%d:%d,%d]" % (p1,p2,row), output=tempfile, format="onedspec")
			os.system("mv %s %s" % (tempfile+".0001.fits", tempfile))
			iraf.noao.twodspec.longslit.scopy(input = err_filename+"[%d:%d,%d]" % (p1,p2,row), output=temp_erfile, format="onedspec")
			os.system("mv %s %s" % (temp_erfile+".0001.fits", temp_erfile))

			# 2. Doppler correct each file.
			redshift_at_aperture = splinez(row+0.5) # 0.5 implies central row
			# Deredshift the spectrum.
			iraf.noao.twodspec.longslit.dopcor(input=tempfile, output=tempfile, redshift = float(redshift_at_aperture))
			iraf.noao.twodspec.longslit.dopcor(input=temp_erfile, output=temp_erfile, redshift = float(redshift_at_aperture))
		
		# 3. Add all the spectra.
		temp_aperture_files = [ "temp%d.fits" % i for i in range(nsum) ]
		temp_aperture_errfiles = [ "temp%d_err.fits" % i for i in range(nsum) ]
		aperture_filename = "Aperture_%d.dat" % aperture 
		aperture_err_filename = "Aperture_%d_err.dat" % aperture 
		SumUpSpectra( temp_aperture_files, aperture_filename, temp_aperture_errfiles, aperture_err_filename)

		# Now, to get yvar from the file.
		yvar = np.loadtxt(aperture_filename, usecols=(1,), unpack=True)
		yvar_err = np.loadtxt(aperture_file, usecols=(2,), unpack=True)

		# Now to do the actual fitting.
		try:
			c, c_err, w, w_err = FitGaussians(xvar, yvar, yvar_err, no_components, rest_wavelengths, mode="both")
			widths.append(w)
			centroids.append(c)
			widths_err.append(w_err)
			centroids_err.append(c_err)
			spatial_points.append(aperture + float(stepsize)/2.0)
		except:
			pass

	for aperture in range(central_row, central_row - uppercut, -1*stepsize):
		plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture ) # Diagnostic
		
		# Not to construct yvar for this aperture. This involves three steps.
		for row in range(aperture, aperture+nsum+1):
			
			# 1. Take each row and throw in a separate fits file.
			tempfile = "temp%d.fits" % counter
			temp_erfile = "temp%d_err.fits" % counter

			iraf.noao.twodspec.longslit.scopy(input = filename+"[%d:%d,%d]" % (p1,p2,row), output=tempfile, format="onedspec")
			os.system("mv %s %s" % (tempfile+".0001.fits", tempfile))
			iraf.noao.twodspec.longslit.scopy(input = err_filename+"[%d:%d,%d]" % (p1,p2,row), output=temp_erfile, format="onedspec")
			os.system("mv %s %s" % (temp_erfile+".0001.fits", temp_erfile))

			# 2. Doppler correct each file.
			redshift_at_aperture = splinez(row+0.5) # 0.5 implies central row
			# Deredshift the spectrum.
			iraf.noao.twodspec.longslit.dopcor(input=tempfile, output=tempfile, redshift = float(redshift_at_aperture))
			iraf.noao.twodspec.longslit.dopcor(input=temp_erfile, output=temp_erfile, redshift = float(redshift_at_aperture))
		
		# 3. Add all the spectra.
		temp_aperture_files = [ "temp%d.fits" % i for i in range(nsum) ]
		temp_aperture_errfiles = [ "temp%d_err.fits" % i for i in range(nsum) ]
		aperture_filename = "Aperture_%d.dat" % aperture 
		aperture_err_filename = "Aperture_%d_err.dat" % aperture 
		SumUpSpectra( temp_aperture_files, aperture_filename, temp_aperture_errfiles, aperture_err_filename)

		# Now, to get yvar from the file.
		yvar = np.loadtxt(aperture_filename, usecols=(1,), unpack=True)
		yvar_err = np.loadtxt(aperture_file, usecols=(2,), unpack=True)

		# Now to do the actual fitting.
		try:
			c, c_err, w, w_err = FitGaussians(xvar, yvar, yvar_err, no_components, rest_wavelengths, mode="both")
			widths.append(w)
			centroids.append(c)
			widths_err.append(w_err)
			centroids_err.append(c_err)
			spatial_points.append(aperture + float(stepsize)/2.0)
		except:
			pass


	widths = np.array(widths)
	centroids = np.array(centroids)
	widths_err = np.array(widths_err)
	centroids_err = np.array(centroids_err)
	spatial_points = np.array(spatial_points)
	filtered = FilterErrorBars(widths_err)

	spatial_points_arc = (spatial_points - centroid)*pixel_scale

	sigma = speed_light*widths/centroids
	sigma_err = sigma*np.sqrt( (centroids_err/centroids)**2 + (widths_err/widths)**2 )

	fig2 = plt.figure(2)
	plt.xlabel("Distance from Center [arcsec]", fontsize=18)
	plt.ylabel("Velocity Dispersion [km/s]", fontsize=18)
	plt.title("For convenience, points with large errors not plotted.")
	for i in range(no_components):
		plt.errorbar( spatial_points_arc[filtered[:,i]], sigma[:,i][filtered[:,i]], sigma_err[:,i][filtered[:,i]], fmt="o" )
		plt.show()
		raw_input()
		plt.clf()

	continue_flag = input_str("Satisfied with at least one rotation curve? (y/n): ")
	if continue_flag == "y":
		# Ask user for file name to store rotation data.
		print("Remember: The outlier points will be saved for completeness, the fitting routine will filter them eventually.")
		filename = raw_input("Enter Base File Name to Save Data: ").rstrip()
		# Output all rotation curve data.
		for i in range(no_components):
			out = Table()
			out.add_column( Column(spatial_points, name="Row") )
			out.add_column( Column(spatial_points_arc, name="DFCG") )
			out.add_column( Column(widths[:,i], name="width") )
			out.add_column( Column(widths_err[:,i], name="width_err") )
			out.add_column( Column(sigma[:,i], name="sigma") )
			out.add_column( Column(sigma_err[:,i], name="sigma_err") )
			out.add_column( Column(sigma_corr[:,i], name="sigmac") )
			out.add_column( Column(sigma_corr_err[:,i], name="sigmac_err") )
			out.write(spectrum_file + filename+"_%d.dat" % (i+1), format="ascii", delimiter="\t")
		
	con = input_str("Do you want another rotation curve? (y/n): ")
	if con == "n":
		break
