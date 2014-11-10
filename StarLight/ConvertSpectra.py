"""
Convert the FITS spectra made by the RotationCurve suite into Starlight
compatible form. 

Usage: python ConvertSpectra.py <2dFitsFile> [clean/backup]

The optional clean statement causes original files to be deleted.
The optional backup statement causes original files to be backed up.
"""

import numpy as np
import pysynphot as ps
from astropy.table import Table
import sys
from SpectralTools import *
import os

try:
	filename = sys.argv[1]
except:
	print("No 2d file name provided.")
	sys.exit(2)

try:
	flag = sys.argv[2]
	if flag == "clean":
		clean=True
	elif flag == "backup":
		backup=True
except:
	clean, backup = False, False

# Load aperture list.
try:
	aperture_table = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
except:
	print("Could not find aperture map. Did you run Derotation algorithm on your file?")
	sys.exit(3)

for spec in aperture_table["col1"]:
	Fits2StarlightAscii(spec)
	if clean:
		os.system("rm %s" % spec)
	if backup:
		try:
			os.system("mkdir %s_backup" % filename[:-5])
			os.system("mv %s %s_backup" % (spec, filename[:-5]) )
		except:
			pass
