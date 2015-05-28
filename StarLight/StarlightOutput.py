"""
This program contains a function to represent the output from Starlight
as a Python object.
"""
from astropy.table import Table
import os
import numpy as np

class ConstructionException(Exception):
	pass

class StarOutput(object):
	"""
	A class to represent the starlight output.
	"""
	def __init__(self, outfilename):
		"""
		Read the output file and initialize some variables.
		"""
		try:
			self.star = open(outfilename, 'r').readlines()
		except:
			print("File needed to construct an output representation object\
					is either absent or corrupt. Please check if file exists.")
			raise ConstructionException

		# Extract all relevant information output by Starlight.
		# 1. No of base spectra and pts in spectra, and goodness-of-fit
		self.no_bases = int(self.star[9].split()[0])
		self.no_points = int(self.star[3*self.no_bases+80].split()[0])
		self.chi2_reduced =  float(self.star[49].split()[0])
		self.adev =  float(self.star[50].split()[0])

		# 2. Kinematic information
		self.vel0 = float(self.star[57].split()[0])
		self.veldisp = float(self.star[58].split()[0])

		# 3. Extinction information
		self.av = float(self.star[59].split()[0])

		# Create an Astropy table based on the above temporary file.
		self.poptable = Table.read(self.star[62:62+self.no_bases], format="ascii")
	
		# Derive some values from the population table.
		# Fractions in Old, Intermediate and Young Populations
		young = (self.poptable["age_j(yr)"] < 1e8)
		inter = (self.poptable["age_j(yr)"] > 1e8) & (self.poptable["age_j(yr)"] < 1e9)
		older = (self.poptable["age_j(yr)"] > 1e9)
		self.old_fraction = np.sum(self.poptable["x_j(%)"][older]) / np.sum(self.poptable["x_j(%)"])
		self.inter_fraction = np.sum(self.poptable["x_j(%)"][inter]) / np.sum(self.poptable["x_j(%)"])
		self.young_fraction = np.sum(self.poptable["x_j(%)"][young]) / np.sum(self.poptable["x_j(%)"])
		
		# Now, read in the observed + model spectrum table.
		self.modelspec = Table.read(self.star[3*self.no_bases+81:3*self.no_bases\
									+81+self.no_points], format="ascii")
		self.modelspec = self.modelspec[ self.modelspec["col4"] != -2 ]
		# Intialization Complete.
	
	def light_meanz(self):
		"""
		Determines the mean metallacity using light vector percentages for weighing.
		"""
		popvector = self.poptable["x_j(%)"]
		metallicities = self.poptable["Z_j"]

		return np.sum(popvector*metallicities) / np.sum(popvector)
	
	def mass_meanz(self):
		"""
		Determines the mean metallacity using mass vector percentages for weighing.
		"""
		popvector = self.poptable["Mini_j(%)"]
		metallicities = self.poptable["Z_j"]

		return np.sum(popvector*metallicities) / np.sum(popvector)

	def light_meanage(self):
		"""
		Determines the mean log age using light vector percentages for weighing. 
		Age output is in log.
		"""
		popvector = self.poptable["x_j(%)"]
		log_ages = np.log10(self.poptable["age_j(yr)"])

		return np.sum(popvector*log_ages) / (np.sum(popvector))

	def mass_meanage(self):
		"""
		Determines the mean log age using light vector percentages for weighing. 
		Age output is in log.
		"""
		popvector = self.poptable["Mini_j(%)"]
		log_ages = np.log10(self.poptable["age_j(yr)"])

		return np.sum(popvector*log_ages) / (np.sum(popvector))

class StarError(object):
	"""
	To represent error information obtained through simulations.
	"""
	def __init__(self, filename):
		er_table = Table.read(filename, format="ascii")

		self.mass_meanage = er_table["StdDev"][9]
		self.light_meanage = er_table["StdDev"][8]
		self.mass_meanz = er_table["StdDev"][7]
		self.light_meanz = er_table["StdDev"][6]
		self.young_fraction = er_table["StdDev"][5]
		self.inter_fraction = er_table["StdDev"][4]
		self.old_fraction = er_table["StdDev"][3]
		self.vel0 = er_table["StdDev"][0]
		self.veldisp = er_table["StdDev"][1]
		self.av = er_table["StdDev"][2]
