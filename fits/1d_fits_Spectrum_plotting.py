#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Viewing a 1d spectrum in fits format,
Read out and display the header data.
Plotting the spectrum.

Version 20230527
@author: lothar schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Please adjust the path and name of the spectrum file
file = input("Enter path and file name, wildcards can be used: ")
file = glob.glob(file)
file = file[0]

#   Reading the spectrum
sp = fits.open(file)

# Read header and print in the console
HD = dict(sp[0].header)
print("\n\nSpectrum header :\n")
for i in HD:
    print(i, ':', HD[i])

try:
    sp[0].header["CRPIX1"]
except:
    sp[0].header["CRPIX1"] = 1

# Creating arrays with the wavelengths and fluxes of the spectrum
flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)
for i in range(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i + 1 - sp[0].header["CRPIX1"]) * sp[0].header["CDELT1"]
    )

#   Closing the fits-file
sp.close()

# Plot full spectrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux, "b-", linewidth=0.5)
plt.xlabel("Wavelength [Angström]", fontsize=14)
plt.ylabel("ADU", fontsize=14)
plt.title("Spektrum " + file, fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.pause(1)

# Showing the plot when the script is executed in a normal Python console
# This is not necessary for the IPython console in Spyder.
# plt.show(block=True)
