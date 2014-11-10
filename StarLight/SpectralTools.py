import numpy as np
import pysynphot as ps
import astropy.io.fits as fits
import astropy.wcs as wcs


def Fits2StarlightAscii(filename):
	"""
	Input: a FITS filename
	Output: None.
	
	The input FITS file is converted into an ASCII form.
	"""
	hdulist = fits.open(filename)
	header = hdulist[0].header
	data = hdulist[0].data
	w = wcs.WCS(header)
	pixels = np.arange(len(data)) + 0.5
	wavelengths = w.all_pix2sky(pixels, 1)[0]

	spectrum = ps.ArraySpectrum(wave=wavelengths, flux=data)
	lowest_wavelength = np.min(spectrum.wave)
	highest_wavelength = np.max(spectrum.wave)

	out_wavelengths = np.arange( np.ceil(lowest_wavelength), np.floor(highest_wavelength), 1)
	flux_out = spectrum.sample(out_wavelengths)
	flux_err = np.sqrt(flux_out)

	f = open(filename[:-5]+".dat","w")
	for i in range(len(flux_out)):
		f.write("  %.0f  %f  %f\n" % (out_wavelengths[i], flux_out[i], flux_err[i]))
	f.close()








