"""
Some support functions for tools such as DerotExtract.py and others.
"""

import pyfits as pf
import numpy as np
from astropy.io import fits
import astropy.wcs as pw
import pysynphot as ps

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


def SarithSpectra( spec_list, output_name, spec_err_list=False, output_err_filename=False, clean=True):
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
	w1 = np.ceil(w1)
	w2 = np.floor(w2)
	wavelengths = np.arange(w1,w2+1,1)

	# Let us start with a loop for spectra.
	spectra_resampled = np.ones( (len(spec_list),len(wavelengths)), float)
	for i,spec in enumerate(spec_list):
		h, d, w = read_1dspectrum(spec)
		s = ps.ArraySpectrum(wave=w, flux=d)
		s = s.resample(wavelengths)
		spectra_resampled[i] = s.flux
	spectrum_mean = spectra_resampled.mean(axis=0)
	
	error_spectra = np.ones( (len(spec_list),len(wavelengths)), float)
	for i,spec_err in enumerate(spec_err_list):
		h, d, w = read_1dspectrum(spec_err)
		s = ps.ArraySpectrum(wave=w, flux=d)
		s = s.resample(wavelengths)
		error_spectra[i] = s.flux**2
	spectrum_mean_error = np.sqrt( error_spectra.mean(axis=0) )

	f = open(output_name, "w")
	for i in range(len(wavelengths)):
		f.write("  %.0f  %f  %f\n" % (wavelengths[i], spectrum_mean[i], spectrum_mean_error[i]))
	f.close()

