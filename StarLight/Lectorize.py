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
	os.system("lector < lector.in")
	os.rename("temp.spec_LINE", spec+"_LINE")
	os.rename("temp.spec_ROSE", spec+"_ROSE")


# Loop over all 1-d spectra.
for spec in aper_map["col1"]:
	# Load spectrum.
	RunLector(spec)
	
# Done.
