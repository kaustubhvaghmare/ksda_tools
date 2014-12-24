"""
A simple diagnostic tool to simply plot the spectrum and the error frame.
"""

import sys
import matplotlib.pyplot as plt
from MatplotlibCustom import *
import numpy as np
from astropy.table import Table

spectrum_file = sys.argv[1]
t = Table.read(spectrum_file, format="ascii")

fig1 = plt.figure(1, figsize=(12,6))
plt.plot(t["col1"], t["col2"], color="blue")
plt.plot(t["col1"], t["col2"]+t["col3"], color="blue", alpha=0.2)
plt.plot(t["col1"], t["col2"]-t["col3"], color="blue", alpha=0.2)
plt.xlabel("Wavelengths [Angstroms", fontsize=18)
plt.ylabel("Relative Flux", fontsize=18)
plt.show()
raw_input()



