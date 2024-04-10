#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import of a series of heliocentrically corrected spectra in fits format,
Display of a selectable line in velocity space in two plots:
    1. all spectra plotted on top of each other.
    2. all spectra with a selectable offset plotted on top of each other.
The plots can be saved as pdf.
The file ListOfSpectralLines.py must be in the same directory as the spectra.

20230216

@author: lothar
"""

import numpy as np
from astropy.io import fits
import glob
import matplotlib.pylab as plt

# local module, must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Line selection
linie, wellenlaenge = Linienlisten.linienauswahl()

bereich = float(input('Enter the speed range to be displayed around \
the selected line in km/s: '))

# Enter the system speed:
systemgeschwindigkeit = float(
    input("Enter the system speed in km/s: "))
systemgeschwindigkeit = systemgeschwindigkeit / 299772 * wellenlaenge

# Parameter for the distance between the spectra in the offset plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")


# Working through the spectra list:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']

    wave = np.zeros(header['NAXIS1'])
    for j in range(header['NAXIS1']):
        wave[j] = lambda0 + j * step

    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

    vstep = step / wellenlaenge * 299772
    pixbereich = int(bereich / vstep)

    fluxbereich = flux[(index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
    wellenlaengenbereich = wave[(
        index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
    geschwindigkeitsbereich = (wellenlaengenbereich - wellenlaenge) / \
        wellenlaenge * 299772

    plt.figure(1, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich, linewidth=.5)
    plt.xlabel('Velocity relative to the rest wavelength in km/s')
    plt.ylabel('relative intensity')
    plt.title(obj + ', ' + linie)
    plt.pause(.01)

    plt.figure(2, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich +
             i * offset, "-", label=JD, linewidth=1)
    # Labeling of the individual spectra:
    plt.text(geschwindigkeitsbereich[1], 1.0 +
             i * offset, format(JD, '.2f'), ha="left", size=7)
    plt.xlabel('Velocity relative to the rest wavelength in km/s')
    plt.ylabel('relative intensity')
    plt.title(obj + ', ' + linie)
    plt.axvline(x=0, color='k', linewidth=.05)
    plt.grid(visible=True, axis='x')
    plt.pause(.01)

frage = input(
    'Would you like to save the graphics as pdf? If yes, enter "y": ')
if frage == 'y':
    plt.figure(1)
    # plt.savefig('Overplot_' + obj + '_'+linie + '.pdf')
    plt.savefig('Overplot_' + obj + '_'+linie + '.png')
    plt.figure(2)
    # plt.savefig('OverplotWithOffset' + obj + '_'+linie + '.pdf')
    plt.savefig('OverplotWithOffset' + obj + '_'+linie + '.png')

plt.close('all')

print('\nEnd of program')
