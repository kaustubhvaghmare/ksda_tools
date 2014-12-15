"""
Usage: python GetCurve.py <2d_spectrum_file> <2d_error_file>
Displays central spectrum.
Allows you to mark wavelengths using m.
Goes and fits Gaussians across the whole range.
"""

import sys
import matplotlib.pyplot as plt
from MatplotlibCustom import *
from RotationCurveLib import *
from helcorr import helcorr
from astropy.table import Table, Column
import pickle

config_file = "GetCurve.config"

# Parse command line arguments.
try:
	spectrum_file = sys.argv[1]
except:
	print("Insufficient number of arguments.")
	print("Usage: python GetCurve.py <2dspec_file> <2d_specErr>")
	sys.exit(2)

try:
	error_file = sys.argv[2]
except:
	error_file = 0

# Check if files can be loaded with a straight face!
header, data, wavelengths = read_spectrum( spectrum_file )
if header == -1:
	print("File corrupt or does not exist.")
	sys.exit(3)
if error_file:
	er_header, er_data, er_wavelengths = read_spectrum( error_file )
	if er_header == -1:
		print("Error file corrupt or does not exist.")
		sys.exit(4)
# DONE! Files now loaded.

dopcor = input_str("Do you want to display spectrum in an approximately deredshifted frame? ")
if dopcor == "y":
	redshift = float(raw_input("Approximate redshift, please? "))
else:
	redshift = 0.0

#############################################
#############################################
# To display the spectrum for spatial centre.
#############################################
#############################################

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
	plt.title("Showing spectrum from row %d of %s" % (central_row, spectrum_file), fontsize=18)
	plt.xlabel("Wavelengths [Angstroms]", fontsize=18)
	plt.ylabel("Relative Flux", fontsize=18)
	plt.grid()
	c = Coordinate(fig1)
	plt.show()
	dumb_input = raw_input()
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
	########################################################
	# Rotation Curve Extraction
	########################################################
	########################################################
	
	# First off, we need to access data from the spectrum but for user
	# convenience we have dealt with wavelengths. Let's get the pixels, shall we?
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
	xvar = np.array(xvar)
	
	# Initialize empty lists to store centroids and error bars in them.
	centroids = [] 
	centroid_errs = []
	spatial_points = []
	for aperture in range(central_row, central_row + uppercut, stepsize):
		plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture )
		yvar = np.mean(data[aperture:aperture+nsum,p1:p2], axis=0)
		if error_file:
			yvar_err = np.mean(er_data[aperture:aperture+nsum,p1:p2], axis=0)
		else:
			yvar_err = None
		try:
			c, cerr = FitGaussians(xvar, yvar, yvar_err, no_components)
			centroids.append(c)
			centroid_errs.append(cerr)
			spatial_points.append( aperture + float(stepsize)/2.0 )
		except:
			pass
	
	for aperture in range(central_row, central_row - uppercut, -1*stepsize):
		plt.plot( np.linspace(x1, x2, p2-p1),  np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) , label="%d" % aperture )
		yvar = np.mean(data[aperture:aperture+nsum,p1:p2], axis=0) 
		if error_file:
			yvar_err = np.mean(er_data[aperture:aperture+nsum,p1:p2], axis=0)
		else:
			yvar_err = None
		try:
			c, cerr = FitGaussians(xvar, yvar, yvar_err, no_components)
			centroids.append(c)
			centroid_errs.append(cerr)
			spatial_points.append( aperture + float(stepsize)/2.0 )
		except:
			pass

	centroids = np.array(centroids)
	centroid_errs = np.array(centroid_errs)
	spatial_points = np.array(spatial_points)

	filtered = FilterErrorBars(centroid_errs)
	filtered_fl = filtered.flatten()
	centroids = centroids[filtered]
	centroids = np.reshape(centroids, (len(centroids),1))
	centroid_errs = centroid_errs[filtered]
	centroid_errs = np.reshape(centroid_errs, (len(centroid_errs),1))
	spatial_points = spatial_points[filtered_fl]

	# Convert spatial row numbers to arc second distance from centre.
	spatial_points_arc = (spatial_points - centroid)*pixel_scale
	
	# Convert wavelength centroids into redshifts.
	#centroids = np.array(centroids)
	#centroid_errs = np.array(centroid_errs)
	centroids_redshifts = (centroids - rest_wavelengths)/rest_wavelengths
	# Compute errors on redshifts.
	centroids_redshifts_err = np.sqrt( (1.0/rest_wavelengths)**2*np.array(centroid_errs)**2 )
	# Convert wavelength centroids into km/s using rest-wavelength comparison. 
	centroids_vel = (np.array(centroids) - rest_wavelengths)/rest_wavelengths * speed_light
	
	# Compute heliocentric correction.
	
	# First, get observatory information.
	lon, lat, alt = ObservatoryInformation(config_file)
	# Then, get object coordinates from input file.
	ra, dec, jd = ObjectInformation(header)
	
	correction, hjd = helcorr(lon, lat, alt, ra, dec, 2455084.537578)
	# Apply heliocentric correction.
	centroids_vel_helio = centroids_vel + correction
	# Transform error bars accordingly.
	centroids_vel_helio_err = np.sqrt( (speed_light/rest_wavelengths)**2*np.array(centroid_errs)**2 + 1e-6)
	
	
	# Plot rotation curve / curves.
	fig2 = plt.figure(2)
	plt.xlabel("Distance from Center [arcsec]", fontsize=18)
	plt.ylabel("Heliocentric Radial Velocity [km/s]", fontsize=18)
	for i in range(no_components):
		plt.errorbar( spatial_points_arc, centroids_vel_helio[:,i], centroids_vel_helio_err[:,i], fmt="o" )
		plt.show()
		raw_input()
		plt.clf()
	
	continue_flag = input_str("Satisfied with at least one rotation curve? (y/n): ")
	if continue_flag == "y":
		# Ask user for file name to store rotation data.
		filename = raw_input("Enter Base File Name to Save Data: ").rstrip()
		# Output all rotation curve data.
		for i in range(no_components):
			out = Table()
			out.add_column( Column(spatial_points, name="Row") )
			out.add_column( Column(spatial_points_arc, name="DFCG") )
			out.add_column( Column(centroids[:,i], name="lambda") )
			out.add_column( Column(centroid_errs[:,i], name="lambda_err") )
			out.add_column( Column(centroids_vel[:,i], name="V") )
			out.add_column( Column(centroids_vel_helio[:,i], name="Vhel") )
			out.add_column( Column(centroids_vel_helio_err[:,i], name="Vhel_rr") )
			out.add_column( Column(centroids_redshifts[:,i], name="z") )
			out.add_column( Column(centroids_redshifts_err[:,i], name="z_err") )
		
			out.write(filename+"_%d.dat" % (i+1), format="ascii", delimiter="\t")
		
		rest_file = open(filename+"_rest.out","w")
		pickle.dump(rest_wavelengths, rest_file)
		rest_file.close()
	
	con = input_str("Do you want another rotation curve? (y/n): ")
	if con == "n":
		break
