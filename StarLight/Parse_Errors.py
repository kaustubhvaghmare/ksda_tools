"""
Parse all the simulations and collect statistics of the various
quantities computed by Starlight.
"""
import glob
import sys
import numpy as np
from StarlightOutput import StarOutput
import os
from astropy.table import Table

# Parse command line arguments.
try:
	filename = sys.argv[1]
except:
	print("No spectrum file specified.")
	sys.exit(2)

# Load aperture map table and get list of aperture spectra names.
aperture_map = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
aper_files = aperture_map["col1"]

# Loop over all apfiles.
for apfile in aper_files:
	# Change to the directory containing simulated outputs.
	os.chdir(apfile[:-4]+'_sim')
	# Gather all realization outputs
	realizations = glob.glob('*.out')

    # Construct arrays to hold the 100 output values.
	vel0s = np.zeros(len(realizations), dtype=float)
	veldisps = np.zeros(len(realizations), dtype=float)
	avs = np.zeros(len(realizations), dtype=float)
	old_fractions = np.zeros(len(realizations), dtype=float)
	inter_fractions = np.zeros(len(realizations), dtype=float)
	young_fractions = np.zeros(len(realizations), dtype=float)
	light_meanzs = np.zeros(len(realizations), dtype=float)
	mass_meanzs = np.zeros(len(realizations), dtype=float)
	light_meanages = np.zeros(len(realizations), dtype=float)
	mass_meanages = np.zeros(len(realizations), dtype=float)
	
	# Loop over all realizations and query the parameters.
	for count, r in enumerate(realizations):
		sl_output = StarOutput(r)
		vel0s[count] = sl_output.vel0
		veldisps[count] = sl_output.veldisp
		avs[count] = sl_output.av
		old_fractions[count] = sl_output.old_fraction
		inter_fractions[count] = sl_output.inter_fraction
		young_fractions[count] = sl_output.young_fraction
		light_meanzs[count]= sl_output.light_meanz()
		mass_meanzs[count] = sl_output.mass_meanz()
		light_meanages[count] = sl_output.light_meanage()
		mass_meanages[count] = sl_output.mass_meanage()
	# Done.
	
	os.chdir('../')
	# Now, compute standard deviation and means, throw them in a file.
	err_output = open(apfile[:-4]+'_er.out', 'w')
	err_output.write("#Quantity Mean StdDev\n")
	err_output.write("Vel0 %.4e %.4e\n" % (np.mean(vel0s), np.std(vel0s)) )
	err_output.write("VelDisp %.4e %.4e\n" % (np.mean(veldisps), np.std(veldisps)) )
	err_output.write("Av %.4e %.4e\n" % (np.mean(avs), np.std(avs)) )
	err_output.write("FracOld %.4e %.4e\n" % (np.mean(old_fractions), np.std(old_fractions)) )
	err_output.write("FracInter %.4e %.4e\n" % (np.mean(inter_fractions), np.std(inter_fractions)) )
	err_output.write("FracYoung %.4e %.4e\n" % (np.mean(young_fractions), np.std(young_fractions)) )
	err_output.write("Lz %.4e %.4e\n" % (np.mean(light_meanzs), np.std(light_meanzs)) )
	err_output.write("Mz %.4e %.4e\n" % (np.mean(mass_meanzs),	np.std(mass_meanzs)) )
	err_output.write("LAge %.4e %.4e\n" % (np.mean(light_meanages), np.std(light_meanages)) )
	err_output.write("MAge %.4e %.4e\n" % (np.mean(mass_meanages), np.std(mass_meanages)) )
