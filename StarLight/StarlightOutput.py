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

		# Now, head towards parsing the output file.
		self.no_bases = int(self.star[9].split()[0])
		self.no_points = int(self.star[3*self.no_bases+80].split()[0])
		self.vel0 = float(self.star[57].split()[0])
		self.veldisp = float(self.star[58].split()[0])

		# Transfer data containing poulations table into a temporary file.
		self.temp_table = open("temp.txt","w")
		self.temp_table.writelines( self.star[62:62+self.no_bases])
		self.temp_table.close()

		# Create an Astropy table based on the above temporary file.
		self.poptable = Table.read("temp.txt", format="ascii")
		# And now delete the file.
		os.remove("temp.txt")

		# Now, read in the observed + model spectrum table.

		# Transfer data containing poulations table into a temporary file.
		self.temp_table = open("temp.txt","w")
		self.temp_table.writelines( self.star[3*self.no_bases+81:3*self.no_bases+81+self.no_points])
		self.temp_table.close()

		# Create an Astropy table based on the above temporary file.
		self.modelspec = Table.read("temp.txt", format="ascii")
		self.modelspec = self.modelspec[ self.modelspec["col4"] != -2 ]
		# And now delete the file.
		os.remove("temp.txt")

		# Intialization Complete.
	
	def mean_metallicity(self):
		"""
		Determines the mean metallacity using population vector percentages for weighing.
		"""
		popvector = self.poptable["x_j(%)"]
		metallicities = self.poptable["Z_j"]

		return np.sum(popvector*metallicities) / np.sum(popvector)

	def mean_age(self):
		"""
		Determines the mean age using population vector percentages for weighing. 
		Age output is in Gyrs.
		"""
		popvector = self.poptable["x_j(%)"]
		ages = self.poptable["age_j(yr)"]

		return np.sum(popvector*ages) / (np.sum(popvector) * 1e9)


		





		
