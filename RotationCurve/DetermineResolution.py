"""
Program to determine resolution by analyzing the arc spectrum of the galaxy.
Inputs: TrimmedArcSpectrum RectifiedGalaxyImage
Output: Width of Gaussian Needed to Fit Spectral Line
"""

import sys
import matplotlib.pyplot as plt
from MatplotlibCustom import *
from RotationCurveLib import *

try:
	arc = sys.argv[1]
	spec = sys.argv[2]
except:
	print("No valid input specified.")
	print("usage: python DetermineResolution.py TrimmedArcSpectrum RectifiedGalaxyImage")

header, data, wavelengths = read_spectrum(arc)
spec_header, spec_data, spec_wavelengths = read_spectrum(spec)
if header == -1 or spec_header == -1:
	print("Files corrupt or do not exist.")
	sys.exit(3)

# Done, files have been loaded. Determine centroid for rectified spectrum.
# Determine the centroid row for displaying a spectrum.
centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )
central_row = int(round(centroid, 0))
redshift=0 # Arc Spectra is always at zero redshift.

# Display Arc Spectrum corresponding to that row.
fig1 = plt.figure(1)
plt.plot( wavelengths/(1+redshift), data[central_row, :], color="blue")
plt.xlim( np.min(wavelengths), np.max(wavelengths) )
plt.title("Showing spectrum from row %d of %s" % (central_row, spec), fontsize=18)
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

rest_wavelengths = np.ones(no_components)
for i in range(no_components):
	rest_wavelengths[i] = float(raw_input("Enter Rest Wavelength for Selected Feature(s): "))
	

# Prepare and do a Gaussian fit.
p1, p2 = Wavelength2Pixels(header, x1, x2)
p1, p2 = np.floor(p1), np.ceil(p2)
pixels = np.arange(p1, p2, 1)
xvar = Pixels2Wavelengths(header, *(pixels+0.5) )
xvar = np.array(xvar)
yvar = data[centroid,p1:p2]
yvar_err = np.ones(len(yvar), dtype=float)
c, cerr = FitGaussians(xvar, yvar, yvar_err, no_components, rest_wavelengths, mode="width")

print("Reported Width = %.3f +- %.3f in wavelength units." % (c,cerr))
