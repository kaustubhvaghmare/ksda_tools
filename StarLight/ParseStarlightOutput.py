"""
Program to determine metallicity gradients, age gradients and any other properties.
It will also produce model vs observed spectrum plots.
"""

import matplotlib.pyplot as plt
import numpy as np
from StarlightOutput import StarOutput
import sys
from astropy.table import Table
from MatplotlibCustom import *
import scipy.stats as ss

def RunningAverage(x,y, bins=6):
	"""
	Input: X and Y and number of bins.
	Output: Xbin_centres, Xerr, Ybin_centres, Yerr
	"""
	mean, binedges, num = ss.binned_statistic(x, y, bins=bins)
	error, binedges, num = ss.binned_statistic(x, y, statistic=np.std, bins=bins)
	xerr = binedges[1] - binedges[0]
	x = (binedges[:-1] + binedges[1:])/2.0
	return x, xerr, mean, error


try:
	filename = sys.argv[1]
except:
	print("No 2d file name provided.")
	sys.exit(2)

# Load aperture list.
try:
	aperture_table = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
except:
	print("Could not find aperture map. Did you run Derotation algorithm on your file?")
	sys.exit(3)

# Temporary, convert pixel distances into arcsec.
aperture_table["col2"] = aperture_table["col2"]*0.254

z = np.zeros(len(aperture_table), dtype=float)
age = np.zeros(len(aperture_table), dtype=float)
vel0 = np.zeros(len(aperture_table), dtype=float)
veldisp = np.zeros(len(aperture_table), dtype=float)

fig1 = plt.figure(1, figsize=(12,6))

for i, spec in enumerate(aperture_table["col1"]):

	filename = spec[:-5]+".out"
	s = StarOutput(filename)
	z[i] =  s.mean_metallicity()
	age[i] = s.mean_age()
	vel0[i] = s.vel0
	veldisp[i] = s.veldisp

	modelspec = s.modelspec

	low, high = np.min(modelspec["col1"]), np.max(modelspec["col1"])

	ax1 = plt.axes([0.1, 0.4, 0.8, 0.5])
	ax1.plot(modelspec["col1"], modelspec["col2"], color="green", label="Observed")
	ax1.plot(modelspec["col1"], modelspec["col3"], color="blue", label="Model")
	ax1.set_ylabel("Relative Flux", fontsize=18)
	ax1.set_xlim(low, high)
	ax1.grid(which="both")
	ax1.set_title("Comparison between observed and modelled spectrum", fontsize=14)

	ax2 = plt.axes([0.1,0.1,0.8,0.3])
	residue = (modelspec["col3"] - modelspec["col2"]) / modelspec["col2"]
	ax2.plot(modelspec["col1"], residue, color="red", label="Residual")
	ax2.set_xlabel("Wavelength in Angstroms", fontsize=18)
	ax2.set_ylabel("Residue", fontsize=18)
	ax2.set_xlim(low, high)
	ax2.grid(which="both")
	fig1.savefig(spec[:-5]+"_mod.png", dpi=150)
	plt.clf()

# Make plots.
fig1 = plt.figure(1)
x, xerr, y, yerr = RunningAverage(np.abs(aperture_table["col2"]), z, bins=len(aperture_table))
plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="bo")
plt.xlabel("Distance from Center [arcsec]", fontsize=18)
plt.ylabel("Metallicity", fontsize=18)
fig1.savefig("z_gradient.png", dpi=150)

fig2 = plt.figure(2)
x, xerr, y, yerr = RunningAverage(np.abs(aperture_table["col2"]), age, bins=len(aperture_table))
plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="bo")
plt.xlabel("Distance from Center [arcsec]", fontsize=18)
plt.ylabel("Age [GYrs]", fontsize=18)
fig2.savefig("age_gradient.png", dpi=150)

fig3 = plt.figure(3)
x, xerr, y, yerr = RunningAverage(np.abs(aperture_table["col2"]), vel0, bins=len(aperture_table))
plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="bo")
plt.xlabel("Distance from Center [arcsec]", fontsize=18)
plt.ylabel("Fitted Shift [kms/sec]", fontsize=18)
fig3.savefig("rc.png", dpi=150)

fig4 = plt.figure(4)
x, xerr, y, yerr = RunningAverage(np.abs(aperture_table["col2"]), veldisp, bins=len(aperture_table))
plt.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="bo")
plt.xlabel("Distance from Center [arcsec]", fontsize=18)
plt.ylabel("Velocity Dispersion [km/s]", fontsize=18)
fig4.savefig("veldisp.png", dpi=150)

