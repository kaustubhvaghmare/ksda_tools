"""
Program which assembles it all together. The final output containing
a lot of output information about the entire galaxy.
"""

import numpy as np
import matplotlib.pyplot as plt
from MatplotlibCustom import *
from astropy.table import Table
from StarlightOutput import StarOutput, StarError
import sys
from matplotlib.backends.backend_pdf import PdfPages

# Get input 2d file from command line arguments.
try:
	filename = sys.argv[1]
except:
	print("No spectrum file specified.")
	sys.exit(2)

# Load aperture map table and get list of aperture spectra names.
aperture_map = Table.read(filename[:-5]+"_aper_map.dat", format="ascii")
aper_files = aperture_map["col1"]
aperture_distances = aperture_map["col3"]
no_apertures = len(aperture_map)

# Initialize blank PDF file for storing output.
models = PdfPages(filename[:-5] + '_models.pdf')

lzs = np.zeros(len(aperture_map), dtype=float)
lages = np.zeros(len(aperture_map), dtype=float)
mzs = np.zeros(len(aperture_map), dtype=float)
mages = np.zeros(len(aperture_map), dtype=float)
lzs_err = np.zeros(len(aperture_map), dtype=float)
lages_err = np.zeros(len(aperture_map), dtype=float)
mzs_err = np.zeros(len(aperture_map), dtype=float)
mages_err = np.zeros(len(aperture_map), dtype=float)

youngfrac = np.zeros(len(aperture_map), dtype=float)
interfrac = np.zeros(len(aperture_map), dtype=float)
oldfrac = np.zeros(len(aperture_map), dtype=float)
youngfrac_err = np.zeros(len(aperture_map), dtype=float)
interfrac_err = np.zeros(len(aperture_map), dtype=float)
oldfrac_err = np.zeros(len(aperture_map), dtype=float)

vel0 = np.zeros(len(aperture_map), dtype=float)
veldisp = np.zeros(len(aperture_map), dtype=float)
avs = np.zeros(len(aperture_map), dtype=float)
vel0_err = np.zeros(len(aperture_map), dtype=float)
veldisp_err = np.zeros(len(aperture_map), dtype=float)
avs_err = np.zeros(len(aperture_map), dtype=float)

fig1 = plt.figure(1, figsize=(8.3,11.7))
# Move over each aperture file.
for i, apfile in enumerate(aper_files):
	# Open the main Starlight output.
	slout = StarOutput(apfile[:-4] + '.out')
	# Open the corresponding error table.
	slerr = StarError(apfile[:-4]+'_er.out')

	lzs[i] = slout.light_meanz()
	lages[i] = slout.light_meanage()
	mzs[i] = slout.mass_meanz()
	mages[i] = slout.mass_meanage()
	lzs_err[i] = slerr.light_meanz
	lages_err[i] = slerr.light_meanage
	mzs_err[i] = slerr.mass_meanz
	mages_err[i] = slerr.mass_meanage

	youngfrac[i] = slout.young_fraction
	interfrac[i] = slout.inter_fraction
	oldfrac[i] = slout.old_fraction
	youngfrac_err[i] = slerr.young_fraction
	interfrac_err = slerr.inter_fraction
	oldfrac_err = slerr.old_fraction
	
	vel0[i] = slout.vel0
	veldisp[i] = slout.veldisp
	avs[i] = slout.av
	vel0_err = slerr.vel0
	veldisp_err = slerr.veldisp
	avs_err = slerr.av
	
	# Start to make of the spectrum, residue and the main plots.
	ax1 = plt.axes( [0.1,0.6,0.8,0.35] )
	ax1.plot(slout.modelspec["col1"], slout.modelspec["col2"], color="green", label='Observed')
	ax1.plot(slout.modelspec["col1"], slout.modelspec["col3"], color="blue", label='Modelled')
	ax1.set_xlabel("Wavelength [Angstroms]", fontsize=16)
	ax1.set_ylabel("Relative Flux", fontsize=16)
	ax1.set_xlim(np.min(slout.modelspec["col1"]), np.max(slout.modelspec["col1"]))
	ax1.legend(loc="upper right", frameon=False)

	ax2 = plt.axes( [0.1,0.3,0.8,0.25] )
	ax2.plot(slout.modelspec["col1"], slout.modelspec["col2"] - slout.modelspec["col3"], color="red")
	ax2.set_xlabel("Wavelength [Angstroms]", fontsize=16)
	ax2.set_ylabel("Residue", fontsize=16)
	ax2.set_xlim(np.min(slout.modelspec["col1"]), np.max(slout.modelspec["col1"]))

	ax3 = plt.axes([0.1,0.25,0.8,0.1],frameon=False)
	#ax3.set_visible(False)
	ax3.xaxis.set_visible(False)
	ax3.yaxis.set_visible(False)
	quantities = [ 
['Dist. From Centre', '%.3f' % aperture_distances[i], '', ''],
['RCS',r'%.3f  ' % (slout.chi2_reduced),  'Avg Dev.',r'%.3f   ' % (slout.adev)],
['LMZ',r'%.3f  %.3f' % (slout.light_meanz(), slerr.light_meanz),  'MMZ',r'%.3f  %.3f' % (slout.mass_meanz(), slerr.mass_meanz)],
['LMAge',r'%.3f  %.3f' % (slout.light_meanage(), slerr.light_meanage),  'MMAge',r'%.3f  %.3f' % (slout.mass_meanage(), slerr.mass_meanage)],
['Vel0',r'%.3f  %.3f' % (slout.vel0, slerr.vel0),  'Veldisp',r'%.3f  %.3f' % (slout.veldisp, slerr.veldisp)],
['Av',r'%.3f  %.3f' % (slout.av, slerr.av),  'YoungFrac',r'%.3f  %.3f' % (slout.young_fraction, slerr.young_fraction)],
['InterFrac',r'%.3f  %.3f' % (slout.inter_fraction, slerr.inter_fraction), 'OldFrac',r'%.3f  %.3f' % (slout.old_fraction, slerr.old_fraction)],
			]
	ax3.table(cellText=quantities)

	
	models.savefig()
	fig1.clf()

models.close()

# Define quantities to be plotted.
quantities = [lzs, mzs, lages, mages, youngfrac, interfrac, oldfrac,
			vel0, veldisp, avs]
quantities_err = [lzs_err, mzs_err, lages_err, mages_err, youngfrac_err, 
			interfrac_err, oldfrac_err, vel0_err, veldisp_err, avs_err]
quant_names = ['Light Weighted Mean Z',
				'Mass Weighted Mean Z',
				'Light Weighted Mean Age',
				'Mass Weighted Mean Age',
				r'Fraction of Young Stars ($<10^8$)',
				r'Fraction of Intermediate Stars ($>10^8, <10^9$)',
				r'Fraction of Old Stars ($>10^{10}$)',
				'Vel0',
				'Dispersion Velocity',
				'Av']

results = PdfPages(filename[:-5] + '_results.pdf')
fig2 = plt.figure(2, figsize=(7.5,6))
for q, q_err, q_name in zip(quantities, quantities_err, quant_names):
	plt.errorbar(aperture_distances, q, q_err, fmt="o", ms=5)
	plt.xlabel("Distance from center [arcsec]", fontsize=16)
	plt.ylabel(q_name, fontsize=16)
	results.savefig()

	fig2.clf()

results.close()

logage_bins = np.linspace(7.5, 10.5, 15+1)
logage_binwidth = logage_bins[1] - logage_bins[0]

histories = PdfPages(filename[:-5] + '_histories.pdf')
fig3 = plt.figure(3)
for i, apfile in enumerate(aper_files):
	# Open the main Starlight output.
	slout = StarOutput(apfile[:-4] + '.out')
	# Get the population table.
	poptable = slout.poptable
	# Obtain and normalize population vector.
	pop_vector_norm = poptable['x_j(%)'] / np.sum(poptable['x_j(%)']) * 100
	mas_vector_norm = poptable['Mini_j(%)'] / np.sum(poptable['Mini_j(%)']) * 100
	# Bin the ages.
	memberships = np.digitize( np.log10( poptable['age_j(yr)'] ), logage_bins )
	rebinned_popvec = np.zeros( 15, dtype=float)
	rebinned_masvec = np.zeros( 15, dtype=float)
	for j in np.arange(1,16):
		rebinned_popvec[j-1] = pop_vector_norm[ (memberships == j) ].sum()
		rebinned_masvec[j-1] = mas_vector_norm[ (memberships == j) ].sum()

	ax1 = plt.axes( [0.1,0.55,0.8,0.40] )
	ax1.bar( logage_bins[:-1], rebinned_popvec, width = logage_binwidth )
	#ax1.set_xlabel('log (Age in Yrs)', fontsize=18)
	ax1.set_title('Distance from Center = %.3f' % aperture_distances[i], fontsize=14)
	ax1.set_ylabel('% (light)', fontsize=18)


	ax2 = plt.axes( [0.1,0.1,0.8,0.40] )
	ax2.bar( logage_bins[:-1], rebinned_masvec, width = logage_binwidth )
	ax2.set_xlabel('log (Age in Yrs)', fontsize=18)
	ax2.set_ylabel('% (mass)', fontsize=18)
	
	histories.savefig()
	fig3.clf()

histories.close()
