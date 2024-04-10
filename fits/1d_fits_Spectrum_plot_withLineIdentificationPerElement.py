#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
View a selectable 1d spectrum in fits format.
Enter the ions/elements whose lines are to be marked in the spectrum.
Graphic can be saved as pdf if desired.

Stand 20230905
author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

# local module, must be in the same directory as the script (or the directory
# of the module in the Python path so that Python can find it)
# The module is located in the blocks folder linelist
import Linienlisten

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Pfad und Name des Spektrumfiles bitte anpassen
file = input("Enter path and file name: ")

#   Reading the spectrum
sp = fits.open(file)

# Read header and print in the console
# HD = dict(sp[0].header)
# print("\n\nHeader des Spektrums :\n")
# print(HD)

if 'CRPIX1' not in sp[0].header:
    sp[0].header["CRPIX1"] = 1

# Creating arrays with the wavelengths and fluxes of the spectrum
flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)
for i in range(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i + 1 - sp[0].header["CRPIX1"]) * sp[0].header["CDELT1"]
    )
# The wave list contains the wavelengths of the pixels
# The corresponding intensities in the flux list

# Close the fits-file
sp.close()

# Plot full spectrum
fig = plt.figure(1, figsize=(14, 10))
plt.plot(wave, flux, "b-", linewidth=0.5)
plt.xlabel("Wavelength [Angström]", fontsize=14)
plt.ylabel("ADU", fontsize=14)
plt.title("Spectrum " + file, fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
plt.pause(.1)

# Selection of the lines or wavelengths
element = Linienlisten.elementauswahl()
for i in element:
    plt.axvline(element[i], color='r')
    plt.text(element[i], 0.3, i)
    plt.pause(.1)

frage = input('If you want to mark more lines, enter y. ')
while frage == 'y':
    element = Linienlisten.elementauswahl()
    for i in element:
        plt.axvline(element[i])
        plt.text(element[i], 0.3, i)
        plt.pause(.1)
    frage = input('If you want to mark more lines, enter y. ')

frage2 = input('Do you want to save the graphic? Then enter y: ')
if frage2 == 'y':
    grafik = file.rsplit('.fit', 1)[0] + '.pdf'
    fig.savefig(grafik)

plt.close('all')
print('Ende of program')
