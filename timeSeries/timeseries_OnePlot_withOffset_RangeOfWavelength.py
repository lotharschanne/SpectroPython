#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a time series of 1d spectra in fits format.
A wavelength range can be selected.
Plots the spectra sections with an offset on top of each other and saves the
the plot as .png and .pdf.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the files (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()


# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")


# Input of the imaged wavelength range
lambda_anfang = float(input('Specify the beginning of the \
wavelength range to be imaged: '))
lambda_ende = float(input('Specify the end of the \
wavelength range to be imaged:'))

# Parameters for the distance between the spectra in the second plot
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

fig = plt.figure(1, (7, 12))
zaehler = 0

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i])
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections:
    step = header["CDELT1"]
    refpix = header["CRPIX1"]

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)

    wave = np.zeros(header["NAXIS1"], dtype=float)
    wave_bereich = np.array([])
    flux_bereich = np.array([])
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step
        if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
            wave_bereich = np.hstack([wave_bereich, wave[k]])
            flux_bereich = np.hstack([flux_bereich, flux[k]])
    plt.plot(wave_bereich, flux_bereich + zaehler * offset, "k-", linewidth=1.)
    plt.text(wave_bereich[1], 1.0 + zaehler *
             offset, filelist[i], ha="left", size=7)
    zaehler = zaehler + 1


# Customize plot properties (grid, labels, etc.):
# plt.ylim(0.2, 1.0 + zaehler * offset + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
# fig.savefig(obj + '_'+ str(lambda_anfang)+'_'+ str(lambda_ende) +
#             '.pdf', format='pdf')
fig.savefig(obj + '_' + str(lambda_anfang)+'_' + str(lambda_ende) +
            '.png', format='png')

plt.pause(.1)

print('Click on the last opened diagram to exit the program.')
plt.waitforbuttonpress(-1)
plt.close('all')
