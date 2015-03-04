"""
A class to read Lector outputs.
We want to supply a 2-d spectrum file.
The class will use the aperture map and create a table
of distance from centers in arc seconds and all measurements.
"""

import numpy as np
from astropy.table import Table
import os
import sys

class MissingApertureMap(Exception):
	pass

class MissingBandFile(Exception):
	pass


class LectorOutput:
	"""
	A class to represent the lector output on the extracted spectra
	from a 2d spectrum in an organized Pythonic way.
	"""
	def __init__(self, spectrum_file):
		# Load aperture map.
		try: 
			print spectrum_file[:-5]+"_aper_map.dat"
			self.ap_table = Table.read(spectrum_file[:-5]+"_aper_map.dat", format="ascii")
		except:
			raise MissingApertureMap

		# Load Band file.
		try:
			self.bands = Table.read("BANDS", format="ascii", data_start=4, header_start=None)
		except:
			raise MissingBandFile

		# Define dfcs
		self.dfcs = self.ap_table["col3"]
		# Define index names
		self.lick_index_names = self.bands["col8"]

		self.indices = {}
		self.indices["dfc"] = self.dfcs
		for i in self.lick_index_names:
			self.indices[i] = []

		for spec in self.ap_table["col1"]:
			for i, lick in enumerate(self.lick_index_names):
				self.indices[lick].append(float(open(spec+"_LINE").read().split()[i+1]))

		for i in self.lick_index_names[1:]:
			self.indices[i] = np.array(self.indices[i], dtype=float)


