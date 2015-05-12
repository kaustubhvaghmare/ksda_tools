"""
This program prepares the Starlight input file for each aperture.
This assumes you have run ConvertSpectra and obtained .dat files
for each aperture's corresponding FITS files. Once done, this program
will prepare the Starlight input file i.e. the grid file which contains
information about which to spectra to fit.
"""

import sys
from astropy.table import Table
import os

try:
	filename = sys.argv[1]
except:
	print("No 2d file name provided.")
	sys.exit(2)

try:
	aperture_table = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
except:
	print("Could not find aperture map. Did you run Derotation algorithm on your file?")
	sys.exit(3)

# Construct some basic parameters needed for making the header of the grid file.
num_fits = len(aperture_table)
current_dir = os.getcwd()

grid_header = """{0:d}                                                [Number of fits to run]
/home/kaustubh/Tools/Starlight/STARLIGHTv04/MilesBase/                 [base_dir]
{1:s}/                                                  [obs_dir]
{1:s}/                                                  [mask_dir]
{1:s}/                                                 [out_dir]
-2007200                                         [your phone number]
4730.0                                           [llow_SN]   lower-lambda of S/N window
4780.0                                           [lupp_SN]   upper-lambda of S/N window
3400.0                                           [Olsyn_ini] lower-lambda for fit
8900.0                                           [Olsyn_fin] upper-lambda for fit
1.0                                              [Odlsyn]    delta-lambda for fit
1.0                                              [fscale_chi2] fudge-factor for chi2
FIT                                              [FIT/FXK] Fit or Fix kinematics
1                                                [IsErrSpecAvailable]  1/0 = Yes/No
0                                                [IsFlagSpecAvailable] 1/0 = Yes/No
""".format(num_fits, current_dir)

# Some basic parameters needed to make the rest of the file.
standard_config = "StCv04.C11.config"
common_mask = filename[:-4]+"mask" #"Masks.EmLines.SDSS.gm"
base_file = "Base.miles.Mun1.30.3Z"
redenning_law = "CCM"
initial_vel = 0
initial_vdisp = 150

# Construct the result
grid_rest = ""
for spec in aperture_table["col1"]:
	grid_rest += "%s   %s   %s   %s   %s   %s   %s   %s\n" % (spec, standard_config,
	base_file, common_mask, redenning_law, initial_vel, initial_vdisp, spec[:-4]+".out")

outfile = open(filename[:-5]+"_grid.in","w")
outfile.write(grid_header + grid_rest)
outfile.close()
