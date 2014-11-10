"""
Program to determine metallicity gradients, age gradients and any other properties.
"""

import matplotlib.pyplot as plt
import numpy as np
from StarlightOutput import StarOutput
import sys
from astropy.table import Table

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

z = np.zeros(len(aperture_table), dtype=float)
age = np.zeros(len(aperture_table), dtype=float)
vel0 = np.zeros(len(aperture_table), dtype=float)
veldisp = np.zeros(len(aperture_table), dtype=float)

for i, spec in enumerate(aperture_table["col1"]):

	filename = spec[:-5]+".out"
	s = StarOutput(filename)
	z[i] =  s.mean_metallicity()
	age[i] = s.mean_age()
	vel0[i] = s.vel0
	veldisp[i] = s.veldisp

# Make plots.
fig1 = plt.figure(1)
plt.scatter(aperture_table["col2"], z)
fig1.savefig("z_gradient.png", dpi=150)

fig2 = plt.figure(2)
plt.scatter(aperture_table["col2"], age)
fig2.savefig("age_gradient.png", dpi=150)

fig3 = plt.figure(3)
plt.scatter(aperture_table["col2"], vel0)
fig3.savefig("rc.png", dpi=150)

fig4 = plt.figure(4)
plt.scatter(aperture_table["col2"], veldisp)
fig4.savefig("veldisp.png", dpi=150)






