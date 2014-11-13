"""
A test program that just takes the output of the DerotExtract.py and plots
the spectra to get a visual feel whether some kind of derotation has happened.
"""

import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.io.fits as fits
import astropy.wcs as wcs
import sys

filename = sys.argv[1]

# Load aperture table.
aperture_table = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")

# Initialize a plot figure.
fig1 = plt.figure(1)

# Loop over all apertures.
for aperture in aperture_table["col1"][:10]:
	# Load spectrum
	spectrum = Table.read(aperture)
	plt.plot(spectrum["WAVELENGTH"], spectrum["FLUX"], label=aperture[:-7])

#legend = plt.legend()
#legend.get_title().set_fontsize('3')
plt.xlim(5750,6000)
plt.ylim(2,10)
plt.grid()
fig1.show()
raw_input()
fig1.savefig("test.png", dpi=150)

	

