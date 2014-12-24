"""
A simple diagnostic tool to simply plot the spectrum and the error frame.
"""

import sys
from astropy.io import fits
import astropy.wcs as pw
import matplotlib.pyplot as plt
from MatplotlibCustom import *
import numpy as np

def read_1dspectrum(filename):
	"""
	Input: Filename
	Output: header, data and wavelengths.
	"""
	hdulist = fits.open(filename)
	wcs = pw.WCS(hdulist[0].header)
	data = hdulist[0].data
	wavelengths = wcs.all_pix2world( np.arange(len(data))+0.5, 1)[0]

	return hdulist[0].header, data, wavelengths

# Get system arguments.
spec1 = sys.argv[1]
spec2 = sys.argv[2]

# Read spectra.
h, d, w = read_1dspectrum(spec1)
he, de, we = read_1dspectrum(spec2)

# Plot spectrum.
fig1 = plt.figure(1, figsize=(12,6))
plt.plot( w, d, color="blue")
plt.plot( w, d+de, color="blue", alpha=0.2)
plt.plot( w, d-de, color="blue", alpha=0.2)
plt.xlabel("Wavelengths [Angstroms", fontsize=18)
plt.ylabel("Relative Flux", fontsize=18)
plt.show()
raw_input()
