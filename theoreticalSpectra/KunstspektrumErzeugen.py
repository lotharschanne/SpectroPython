#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Skript erzeugt ein künstliches, per gaußbroadening auf eine bestimmte FWHM,
ausgedrückt durch ein R, verbreitertes Linien-Spektrum

Created on Fri Feb 18 15:38:26 2022

@author: lothar
"""

from PyAstronomy import pyasl
import matplotlib.pylab as plt
import numpy as np
from astropy.io import ascii, fits

# Set up an input spectrum
x = np.linspace(5000.0, 5100.0, 1000)
y = np.ones(x.size)

# Introduce some delta-peaked lines
y[165] = 0.7
y[187] = 0.3
y[505] = 0.1
y[610] = 0.1
y[615] = 0.7
y[900] = 0.1

print("Linien bei: ", x[165], x[187], x[505], "\n", x[610], x[615], x[900])
X = [x[165], x[187], x[505], x[610], x[615], x[900]]

# Apply Gaussian instrumental broadening, setting the resolution to 10000.
r, fwhm = pyasl.instrBroadGaussFast(x, y, 10000, edgeHandling="firstlast", fullout=True)

print("FWHM used for the Gaussian kernel: ", fwhm, " A")

# Plot the output
plt.step(x, r, "r-", label="Broadened curve (full)")
plt.step(x, y, "b-", label="Input")
plt.legend(loc=4)
plt.show()

ascii.write(
    [x, r], "Kunstspektrum.dat", overwrite=True, names=["WAVE", "FLUX"], format="tab"
)

f = open("Kunstspektrum_Linien", "w")
for i in X:
    f.write(str(i) + "\n")
f.write("FWHM: " + str(fwhm))
f.close()

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRVAL1"] = 5000.0
header["NAXIS1"] = 1000
header["CRPIX1"] = 1.0
header["CDELT1"] = 0.1
fits.writeto("Kunstspektrum.fits", r, header, overwrite=True, output_verify="silentfix")
