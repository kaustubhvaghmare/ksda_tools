"""
Program to run Lector on all the extracted spectra from DerotExtract.py.
The program will do the following:
1) Translate the spectra into Lector compatible forms.
2) Run Lector on them.
3) Make a unified table of outputs
"""

import sys
from astropy.table import Table
import os
import numpy as np
import numpy.random as nr

# Parse command line arguments.
try:
	spec_file = sys.argv[1]
except:
	print("Incorrect usage. Specfile should be provided.")
	sys.exit(2)

# Open the aperture map file.
try:
	aper_map = Table.read(spec_file[:-5]+"_aper_map.dat", format="ascii")
except:
	print("You have not run DerotExtract utility on this file. Please do so before using this tool.")
	sys.exit(3)


def RunLector(specfile):
	"""
	Input: Spectrum Astropy Table
	Output: None, but the spectrum is converted into a Lector compatible form and Lector is executed.
	"""
	# Make lector configuration file.
	lect_conf = open("lector.in", "w")
	lect_conf.write("s\ntemp.spec\n0\nn")
	lect_conf.close()

	# Convert spectrum into Lector form.
	spectrum = Table.read(specfile, format="ascii")
	lec_spec = spectrum["col1", "col2"]
	lec_spec.write("temp.spec", format="ascii")

	# Run Lector.
	os.system("/home/kaustubh/Tools/Lector/lector < lector.in > /dev/null")
	os.rename("temp.spec_LINE", spec+"_LINE")
	os.rename("temp.spec_ROSE", spec+"_ROSE")

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
	
	print("Simulating Synthetic Spectra to Obtain Errors")
	for i in range(sim):

		# Make realization.
		noise = nr.normal(0,1,len(spectrum)) * spectrum["col3"]
		realization = spectrum["col1", "col2"]
		realization["col2"] = realization["col2"] + noise
		realization.write("temp.spec", format="ascii")

		# Run lector.
		os.system("/home/kaustubh/Tools/Lector/lector < lector.in > /dev/null")

		sim_data[i] = np.array(open("temp.spec_LINE").read().split()[1:-3], dtype=float)
	
	errors = np.std(sim_data, axis=0)
	err_out = open(spec+"_LINE_ERR", "w")
	err_out.write("temp.spec    ")
	err_out.write("  ".join([str(i) for i in errors]))
	err_out.close()
	os.system("rm temp.spec*")


# Loop over all 1-d spectra.
for spec in aper_map["col1"]:
	print("Working on %s" % spec)
	RunLector(spec) # Runs Lector
	LectorError(spec) # Determines error
	
# Done.
