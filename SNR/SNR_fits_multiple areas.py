#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loads a 1d spectrum in fits format and plots it.
Selection of spectrum sections by 2 mouse clicks and then
pressing the escape button twice.
The SNR of the areas are printed out and saved in a comma-separated csv file
file.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend


# Path and name of the fits-file
file = input("Path and name of the fits-file: ")

#   Read header and data (flux)
flux, header = fits.getdata(file, header=True)

print("Minimum and Maximum in flux: ", flux.min(), "  ", flux.max())


#   Check for the necessary header entries
nax = header["NAXIS1"]
crval = header["CRVAL1"]
cdel = header["CDELT1"]
if "CRPIX1" not in header:
    header["CRPIX1"] = 1

#   Generation of a numpy array with the wavelengths of the spectrum
wave = np.ones(nax, dtype=float)
crval = crval + (1 - header["CRPIX1"]) * cdel
for i in range(nax):
    wave[i] = crval + i * cdel
# The wave list contains the wavelengths of the pixels.
# In the list flux the corresponding intensities.

# Plot spectrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angstroem]")
plt.ylabel("ADU")
plt.title("Spectrum " + file)
plt.grid(True)
plt.pause(10)  # Time to enlarge the desired spectrum section

# Interactive setting of the base points begin and end for SNR estimation
# Press the ESC key twice to end the interactivity
# plt.waitforbuttonpress()

snr = []
Anfang = []
Ende = []
eingabe = "y"
while eingabe == "y":
    print("After marking the two wavelengths in the diagram using the left ")
    print("Press the escape key twice: ")
    pts = []
    pts = np.asarray(plt.ginput(n=-1, timeout=-1))
    if plt.waitforbuttonpress():
        pass
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=6)

    a = pts[0, 0]
    b = pts[1, 0]

    aindex = int((a - crval) / cdel)
    bindex = int((b - crval) / cdel)

    newflux = flux[aindex:bindex]
    newwave = wave[aindex:bindex]

    # Calculation of SNR:
    SNR = newflux.mean() / newflux.std()
    snr.append(SNR)
    Anfang.append(a)
    Ende.append(b)
    print("SNR between %.2f und %.2f  = %.1f" % (a, b, SNR))
    print("\nNEnter a new range? Enter y for this.")
    print("Enter a different letter to continue the program: ")
    eingabe = input()
    if eingabe != "y":
        break

# Saving as an ascii file
ascii.write(
    [Anfang, Ende, snr],
    file + "_SNR.csv",
    overwrite=True,
    names=["Begin", "End", "SNR"],
    format="csv",
)

SNR = 0
for i in range(len(snr)):
    SNR += snr[i]
SNR_mean = SNR / len(snr)
print("Mean SNR = ", SNR_mean)

plt.close('all')
