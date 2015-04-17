"""
Program to use a GUI interface to remove outlier points in a curve or replace
them with interpolations between the two.
"""

from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.interpolate import interp1d
from RotationCurveLib import *

class Coordinate:
	def __init__(self, figure):
		self.coord = []
		self.figure = figure
		figure.canvas.mpl_connect("key_press_event", self)

	def __call__(self, event):
		if event.key == "i":
			self.coord.append(event.xdata)
			print "Coordinate added = %.2f" % (event.xdata)
		elif event.key == "q":
			plt.close(self.figure)

	def getcoords(self):
		return self.coord

def DeletePoint(x, curve):
	"""
	Input: the DFCG point to be deleted and curve from where it has to be deleted.
	Output: Modified curve
	Notes: Curve is assumed to be astropy table object.
	"""
	min_arg = np.argmin(np.abs(curve["DFCG"] - x))
	curve.remove_row(min_arg)
	return curve

def DeletePoints(x, curve):
	for point in x:
		curve = DeletePoint(point, curve)
	return curve

#############################################
#############################################

def InterpFill(x, curve):
	"""
	Input: the DFCG point to be reinterpolated and curve from where it resides!
	Output: Modified curve
	Notes: Curve is assumed to be astropy table object.
	"""	
	min_arg = np.argmin(np.abs(curve["DFCG"] - x))
	target = curve["DFCG"][min_arg]
	columns = ["lambda", "V", "Vhel", "z"]
	for c in columns:
		f = interp1d(curve["DFCG"][[min_arg-1,min_arg+1]], curve[c][[min_arg-1,min_arg+1]], kind="linear")
		curve[c][min_arg] = f(curve["DFCG"][min_arg])
	
	return curve

def InterpFills(x, curve):
	for point in x:
		curve = InterpFill(point, curve)
	return curve
		
try:
	curve_file = sys.argv[1]
	curve = Table.read(curve_file, format="ascii")
	#filtered1 = FilterErrorBars(curve["Vhel_err"])
	#filtered2 = MeanSigmaClipper(curve["Vhel"], sigma=2.50)
	#filtered = (filtered1 & filtered2)
	curve_temp = curve#[filtered]

except:
	print("No or bad curve specified. Try again.")
	sys.exit(1)

fig1 = plt.figure(1)
plt.errorbar( curve_temp["DFCG"], curve_temp["Vhel"], curve_temp["Vhel_err"], fmt="o" )
plt.xlabel(" Distance from Center", fontsize=18)
plt.ylabel("Velocity", fontsize=18)
plt.title("Select point using 'i'", fontsize=18)
plt.grid()
c=Coordinate(fig1)
plt.show()
dumb = raw_input()
x = c.getcoords()

choice = int(raw_input(("1) Delete\n2) InterpFill\nChoice:")))
if choice == 1:
	curve = DeletePoints(x, curve)
elif choice == 2:
	curve = InterpFills(x, curve)
else:
	print("InvalidResponse")
	sys.exit(3)

curve.write(curve_file+"_mod.dat", format="ascii")



