# Program accepts a rotation curve file generated by the GetCurve.py routine
# And the fits a high order spline function to the same.

from astropy.table import Table
import matplotlib.pyplot as plt
from MatplotlibCustom import *
import numpy as np
import sys
import scipy.interpolate as si
import pickle
from numpy.polynomial.polynomial import polyfit, polyval

# Check if file name was supplied.
try:
	filename = sys.argv[1]
except:
	print "Improper syntax. No file name specified."
	sys.exit(2)
# Attempt loading the file.
try:
	rotcurve = Table.read(filename, format="ascii")
except:
	print "No such curve exists or curve corrupt."
	sys.exit(3)

# Display the rotation curve for user to preview.
fig1 = plt.figure(1)
plt.errorbar( rotcurve["Row"], rotcurve["lambda"], rotcurve["lambda_err"], fmt="o" )
plt.xlabel("Row Number", fontsize=18)
plt.ylabel("Wavelength Center", fontsize=18)
plt.title("A spline fit will be performed to the following curve", fontsize=14)
fig1.show()
dum = raw_input("Press any key to continue fitting.")

# A list of types.
#interp_types = ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']


# Get sorted element indices for x-axis variable. It is the requirement of the
# function being used for Spline Fitting that the X be in an increasing order.
# But the way the rotation curve was extracted was to stepping away from the centre
# to the outside and then back to centre and the other side. Hence the need for
# sorting the same.

sortind = np.argsort( rotcurve["Row"] )

fig1 = plt.figure(1)
xs = np.linspace( np.min(rotcurve["Row"]), np.max(rotcurve["Row"]), 1000)
for i in range(10):
	interp_coeff = polyfit( rotcurve["Row"], rotcurve["lambda"], w=1/rotcurve["lambda_err"], deg=i)
	interp_coeff_z = polyfit( rotcurve["Row"], rotcurve["z"], w=1/rotcurve["z_err"], deg=i)


	# Below, we are just generating a sample curve for the user to visualize the fit.
	ys = polyval(xs, interp_coeff)

	# Now, present the plot to the user.
	plt.errorbar( rotcurve["Row"], rotcurve["lambda"], rotcurve["lambda_err"], fmt="o" )
	plt.plot(xs, ys)
	plt.xlabel("Row Number", fontsize=18)
	plt.ylabel("Wavelength Center", fontsize=18)
	plt.title("Fit for degree %d is Shown." % i, fontsize=14)
	fig1.show()
	
	# Formally confirm that the user is satisfied with the spectrum
	choice = raw_input("Are you okay with the fit? 'y' to approve, any other key to try another fit.").rstrip()
	if choice == "y":
		outfile = open( filename+"_fit.out", "w")
		outfilez = open( filename+"_fit_z.out", "w")
		pickle.dump(spline, outfile)
		pickle.dump(spline_z, outfilez)
		print "Best-fit spline saved as %s." % (filename+"_spline.out")
		print "Note that this file can only be retrieved by the derotation program."
		print "Using outside of Derotation program is discouraged."
	plt.clf()
