"""
Program to compute errors on indices returned by Lector.
This involves simulations of the spectra using error information
Rerunning Lector over and over again.
"""

import sys
from astropy.table import Table
import numpy.random as nr
import os
import numpy as np
import numpy.random as nr


# Define a function to determine error bars.
def LectorError(specfile, sim=100):
	"""
	Input: Spectrum Astropy Table
	Output: Errors on indices
	"""
	# Make lector configuration file.
	lect_conf = open("lector.in", "w")
	lect_conf.write("s\ntemp.spec\n0\nn")
	lect_conf.close()
	
	no_bands = len(Table.read("BANDS", format="ascii", data_start=4, header_start=None))
	sim_data = np.zeros( (sim, no_bands), dtype=float)

	# Load the spectrum.
	spectrum = Table.read(specfile, format="ascii")
#	lec_spec = spectrum["col1", "col2"]
	
	for i in range(sim):

		# Make realization.
		noise = nr.normal(0,1,len(spectrum)) * spectrum["col3"]
		realization = spectrum["col1", "col2"]
		realization["col2"] = realization["col2"] + noise
		realization.write("temp.spec", format="ascii")

		# Run lector.
		os.system("/home/kaustubh/Tools/Lector/lector < lector.in > /dev/null")

		sim_data[i] = np.array(open("temp.spec_LINE").read().split()[1:-3], dtype=float)
	
	return np.std(sim_data, axis=0)



