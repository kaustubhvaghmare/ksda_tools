"""
Creates 100 synthetic spectra for all Nap aperture spectra of
the galaxies.
"""

import numpy.random as nr
from astropy.table import Table
import sys
import os

nsim = 100

# Parse command line arguments.
try:
	filename = sys.argv[1]
except:
	print("No spectrum file specified.")
	sys.exit(2)

grid_template = """{0:d}                                                [Number of fits to run]
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
"""

standard_config = "StCv04.C11.config"
common_mask = filename[:-4]+"mask" #"Masks.EmLines.SDSS.gm"
base_file = "Base.miles.Mun1.30.3Z"
redenning_law = "CCM"
initial_vel = 0
initial_vdisp = 150

def WriteTable(table, filename):
	"""
	Input: Astopy Table object representing spectrum and name of output file.
	Output: File compatible with Starlight Spectrum Input.
	"""
	f = open(filename, "w")
	for i in range(len(table)):
		f.write("  %.0f  %f  %f\n" % (table["col1"][i], table["col2"][i], table["col3"][i]) )

# Load the aperture map table.
aperture_map = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
aper_files = aperture_map["col1"]

for apfile in aper_files:
	actual_spectrum = Table.read( apfile, format="ascii")
	os.system("mkdir %s_sim" % apfile[:-4])
	for i in range(nsim):
		noise = nr.normal(0,1,len(actual_spectrum)) * actual_spectrum["col3"]
		realization = actual_spectrum["col2"] + noise
		sim = actual_spectrum.copy()
		sim["col2"] = realization
		WriteTable(sim, "%s_sim/realization_%d.dat"% (apfile[:-4], i+1))
		#sim.write("%s_sim/realization_%d.dat" % (apfile[:-4], i+1), format="ascii.no_header")
	os.chdir("%s_sim" % apfile[:-4])
	star_grid = open("%s_simgrid.in" % apfile[:-4], "w")
	current_dir = os.getcwd()
	star_grid.write( grid_template.format(nsim,current_dir) )
	grid_rest = ""
	for i in range(nsim):
		grid_rest += "%s   %s   %s   %s   %s   %s   %s   %s\n" % ("realization_%d.dat" % (i+1), standard_config,
		base_file, common_mask, redenning_law, initial_vel, initial_vdisp, "realization_%d.out" % (i+1))
	star_grid.write(grid_rest)
	os.system("cp ../StCv04.C11.config ../%s ../Base.miles.Mun1.30.3Z . " % (filename[:-4]+"mask" ) )
	os.system("/home/kaustubh/Tools/Starlight/STARLIGHTv04/StarlightChains_v04.exe < %s_simgrid.in" % apfile[:-4])
	break


		

