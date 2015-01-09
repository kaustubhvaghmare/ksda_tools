"""
Program to make mean rotation curve based on rotation curves obtained for
various features.
"""

import sys
import matplotlib.pyplot as plt
from MatplotlibCustom import *
import glob
from astropy.table import Table, Column
from RotationCurveLib import *

# Parse command line arguments, just one needed.
try: 
	spectrum_file = sys.argv[1]
except:
	print("No 2d spectrum file provided.")
	sys.exit(2)

# Search for rotation curve files.
curve_files_all = glob.glob(spectrum_file+"*.dat")
if len(curve_files_all) == 0:
	print("No rotation curves found for this galaxy. Please run GetCurve first and save some curves.")
	sys.exit(3)

# Load them as table objects.
curve_tables_all = []
for curve in curve_files_all:
	curve_tables_all.append( Table.read(curve, format="ascii") )

# Display the curves one by one.
fig1 = plt.figure(1, figsize=(6,6))
for i,t in enumerate(curve_tables_all):
	filtered = FilterErrorBars(t["z_err"])
	plt.errorbar( t["DFCG"][filtered], t["z"][filtered], t["z_err"][filtered], fmt="o" )
	plt.xlabel("Row Number", fontsize=18)
	plt.ylabel("Wavelength Center", fontsize=18)
	plt.title("For file %s" % curve_files_all[i])
	plt.show()
	raw_input("Press any key for next curve.")
	plt.clf()

print("Curves have not been displayed. You are now required to reject the curves which you believe" 
		"have not been determined sufficiently well. Only the remaining curves, if any, will be used for"
		"constructing the mean curve.")
for i,f in enumerate(curve_files_all):
	print("%d : %s" % (i+1,f))

while True:
	rejected = raw_input("Enter a comma separated list of numbers of files to be rejected: ")
	try:
		rejected = [ int(i)-1 for i in rejected.split(",") ]
		break
	except:
		print("Invalid entries.")
		pass

curve_files = []
curve_tables = []

# Remove rejected
for i in rejected:
	curve_files.append( curve_files_all[i] )
	curve_tables.append( curve_tables_all[i] )

# Next, make the mean curve.
columns = ["lambda", "V", "Vhel", "z"]
err_columns = ["lambda_err", "Vhel_err", "z_err"]
mean_table = Table()
mean_table.add_column( curve_tables[0]["Row"] )
mean_table.add_column( curve_tables[0]["DFCG"] )

for c in columns:
	temp = curve_tables[0][c]
	for i in range(1,len(curve_tables)):
		temp += curve_tables[i][c]
	temp = temp / len(curve_tables)
	mean_table.add_column( temp )

for c in err_columns:
	temp = curve_tables[0][c]**2
	for i in range(1,len(curve_tables)):
		temp += curve_tables[i][c]**2
	temp = temp / len(curve_tables)
	temp = np.sqrt(temp)
	col = Column( temp, name=c)
	mean_table.add_column( col )

# Mean table has been constructed.
# Now, save mean_table.
# Perform outlier rejection on mean curve.
filtered = FilterErrorBars(mean_table["z_err"])
# Filter the mean_table.
mean_table = mean_table[filtered]

# Mean curve has been obtained.
fig1 = plt.figure(1, figsize=(6,6))
plt.errorbar( mean_table["DFCG"], mean_table["z"], mean_table["z_err"], fmt="o" )
plt.xlabel("Row Number", fontsize=18)
plt.ylabel("Wavelength Center", fontsize=18)
plt.title("For file %s" % curve_files[i])
plt.show()
raw_input("Press a key to continue.")
plt.clf()

happy = input_str("Happy with the mean curve? (y/n): ")
if happy == "y":
	mean_table.write(spectrum_file+"_mean.dat", format="ascii")
	print("Curve outputted to file: %s" % (spectrum_file+"_mean.dat"))
else:
	print("You may want to rerun GetCurve or MakeMeanCurve again.")
