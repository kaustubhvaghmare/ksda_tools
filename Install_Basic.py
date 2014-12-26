"""
Simple script that attempts to set up the necessary stuff to make
the basic SALT reduction tools available on a machine.
"""

import sys
import os
import shutil as sh
import glob

# Start by checking installations of third party libraries.
missing = []
try:
	import numpy
except:
	missing.append("numpy")

try:
	import astropy
except:
	missing.append("astropy")

try:
	import pyfits
except:
	missing.append("pyfits")

try:
	import matplotlib
except:
	missing.append("matplotlib")

if missing:
	print("Following packages are needed for tools to work, not found.")
	for i in missing:
		print(i)
	sys.exit(2)
else:
	print("All dependencies satisfied.")

# If program reaches here, it means that dependencies have been satisfied
# and so we can start the process of installation.
home = os.getenv('HOME')
# Query installation directory.
install_dir = (raw_input("Please enter installation path: ")).strip()
if install_dir[-1] == "/": install_dir = install_dir[:-1]

# Verify if directory exists.
if os.path.isdir(install_dir):
	delete = raw_input("Directory already exists. Replace it? (y/n): ").lower()
	if delete == "y":
		# Deleting
		try:
			sh.rmtree(install_dir)
		except:
			print("Could not access the directory. Please check permissions.")

# Copy files to the destination.
print("Copying files to the destination...")
cwd = os.getcwd()
sh.copytree(cwd+"/BasicReduction", install_dir+"/")
print("Done.")

# Make modules executable.
modules = ['Background','Combine','Extract','Flat','Fluxcal','Foreground',\
			'Identify','Preprocess','Refresh','Transform','Understand']
print("Converting scripts to executable form...")
for m in modules:
	f = install_dir+"/"
	os.chmod(f,0744)
print("Done.")

# Home IRAF directory.
print("We need to make changes to your $HOME/iraf directory.")
# Create backup of current.
if os.path.exists(home+"/iraf"):
	sh.copytree(home+"/iraf", home+"/iraf.bak")
	sh.rmtree(home+"/iraf")
	print("Your current IRAF directory backed up into %s" % (home+"/iraf.bak"))
	sh.copytree(install_dir+"/iraf", home)
else:
	sh.copytree(install_dir+"/iraf", home+"/iraf")

# Initialize login.cl for new directory.
print("Running mkiraf.")
os.chdir(home+"/iraf")
os.system("mkiraf")

# Offer .bashrc changes.
print("It's necessary for the installation directory to be in your PATH.")
bash = raw_input("I can modify the .bashrc to this effect. Should I (y/n): ").lower()
if bash == "y":
	with open(home+'/.bashrc','a') as b:
		entry = "\nexport PATH=$PATH:%s" % (install_dir)
		b.write(entry)
	print("A path entry has been updated to bashrc. You may need to source it or restart the shell.")
else:
	print("You need to change your PATH variable to point to the installation directory. Else tools will not work!")

# Finally, provide the irafhome.
with open(install_dir+"/IRAFHome.py","w") as f:
	f.write("irafhome = '%s'" % (install_dir+"/iraf"))
