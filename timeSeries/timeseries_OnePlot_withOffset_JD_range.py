#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a time series of 1d spectra in fits format. A JD range to be plotted
is selected.
Plots the spectra with an offset on top of each other and saves the plot as a
.png and .pdf.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob


# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the file (use wildcards) : ")

filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

JD = np.zeros(len(filelist))
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD[i] = header['JD']
print('JD range from', JD.min(), ' bis ', JD.max())

# Entering the range of observation times as JD
JD_anfang = float(input('Enter the start of the range of \
observation times (JD): '))
JD_ende = float(input('Enter the end of the range of \
observation times (JD): '))

# Parameters for the distance between the spectra in the second plot
offset = float(input("Please enter the desired offset: "))
# obj = input("Please enter the object name: ")

obj = 'Algol JD' + str(JD_anfang) + ' bis ' + str(JD_ende)

fig = plt.figure(1, figsize=(7, 10))
xlim_a = 4330  # Adjust !!!!!!
xlim_e = 4350  # adjust !!!!!!
plt.xlim(xlim_a, xlim_e)
zaehler = 0
# Read header and flux
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD = header['JD']
    JD_float = float(JD)
    if JD_float >= JD_anfang and JD_float <= JD_ende:
        print(filelist[i], ":")
        step = header["CDELT1"]
        refpix = header["CRPIX1"]
        wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)
        zusatz = max(flux) - 1
        wave = np.zeros(header["NAXIS1"], dtype=float)
        for k in range(header["NAXIS1"]):
            wave[k] = wave_erstesPix + k * step
        plt.plot(wave, flux + zaehler * offset, "-", linewidth=1)
        # Labeling of the individual spectra:
        plt.text(xlim_a, 1 + zaehler * offset, JD, ha="left", size=8)
        # plt.xlim(6300,6500)
        plt.pause(0.1)
        zaehler += 1
    else:
        pass


# Customize plot properties (grid, labels, etc.):
# plt.ylim(0.2, 1.0 + (zaehler * offset) + zusatz)
plt.title("Time Series of " + obj)
plt.grid(True)
plt.xlabel("Wavelength in Angström")
plt.ylabel('normierter Flux')
plt.pause(0.1)


frage = input(
    'Would you like to save the graphic? If yes, enter "y":')
if frage == 'y':
    plt.savefig('ofsettedOverplot_' + obj + '.pdf')
    plt.savefig('ofsettedOverplot_' + obj + '.png')

e = input('Press the e button to end the program.')
if e == 'e':
    plt.close('all')
    print('End of the program')
