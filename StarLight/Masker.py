"""
Program to use a GUI interface to mark regions to be masked.
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d
from RotationCurveLib import *

# Parse command line arguments.
try:
	spectrum_file = sys.argv[1]
except:
	print("No spectrum provided.")
	sys.exit(2)

# Read Spectrum
header, data, wavelengths = read_spectrum( spectrum_file )

# Centroid
centroid =  GetCentroid( range(header["NAXIS2"]), np.sum(data, axis=1) )
central_row = int(round(centroid, 0))

class Masker:
	def __init__(self, figure):
		self.counter = 0
		self.coord = []
		self.figure = figure
		figure.canvas.mpl_connect("key_press_event", self)

	def __call__(self, event):
		if event.key == "m":
			if self.counter %2 == 0:
				self.temp = [0,0]
				self.temp[0] = event.xdata
				print "Boundary edge added. Add second edge.",
				self.counter += 1
			else:
				print "..Second edge added."
				self.temp[1] = event.xdata
				self.coord.append(tuple(self.temp))
				ax = self.figure.gca()
				yl, yu = ax.get_ylim()
				meany = (yu+yl)/2.0
				ax.plot( self.temp, [meany,meany], color="red", linewidth=5, alpha=0.5 )
				self.counter+=1
		if event.key == "q":
			if self.counter%2 != 0:
				print("You can't quit now, a boundary edge needs to be added!")
			else:
				plt.close(self.figure)

	def getcoords(self):
		return self.coord

def MakeMask(spectrum_file, regions):
	mask_out = open(spectrum_file[:-4]+"mask", "w")
	mask_out.write("%d\n" % len(regions))
	for reg in regions:
		left = round(np.min(reg))
		right = round(np.max(reg))
		mask_out.write("%.1f  %.1f  0.0  USER MASK\n" % (left, right))
	mask_out.close()


#####################################
# MAIN PROGRAM
#####################################

dopcor = input_str("Do you want to display spectrum in an approximately deredshifted frame? ")
if dopcor == "y":
	redshift = float(raw_input("Approximate redshift, please? "))
else:
	redshift = 0.0

fig1 = plt.figure(1)
plt.plot( wavelengths/(1+redshift), data[central_row, :], color="blue")
plt.xlim( np.min(wavelengths), np.max(wavelengths) )
plt.title("Showing spectrum from row %d of %s" % (central_row, spectrum_file), fontsize=18)
plt.xlabel("Wavelengths [Angstroms]", fontsize=18)
plt.ylabel("Relative Flux", fontsize=18)
plt.grid()
c = Masker(fig1)
plt.show()
dumb_input = raw_input()

regions = c.getcoords()

# Offer a preview of the regions...
fig2 = plt.figure(2)
plt.plot( wavelengths/(1+redshift), data[central_row, :], color="blue")
for reg in regions:
	p1, p2 = Wavelength2Pixels(header, reg[0], reg[1])
	p1,p2 = np.floor(p1), np.ceil(p2)
	y = np.ones(2)*np.mean(data[central_row, p1:p2])
	plt.plot( reg, y , color="red", linewidth=5, alpha=0.5)
plt.xlim( np.min(wavelengths), np.max(wavelengths) )
plt.title("Showing spectrum from row %d of %s" % (central_row, spectrum_file), fontsize=18)
plt.xlabel("Wavelengths [Angstroms]", fontsize=18)
plt.ylabel("Relative Flux", fontsize=18)
plt.grid()
plt.show()
accept = raw_input("Finalize mask? (y/n): ").lower()

if accept == "y":
	MakeMask(spectrum_file, regions)


