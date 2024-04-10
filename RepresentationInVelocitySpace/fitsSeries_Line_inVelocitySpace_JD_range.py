#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in a series of heliocentrically corrected spectra in fits format,
Input of a time period as JD start and JD end.
Display of a selectable line in velocity space in two plots:
    First, all spectra plotted on top of each other.
    Secondly, all spectra with a selectable offset plotted on top of each other.
The plots can be saved as pdf and png.

The file ListOfSpectralLines.py must be in the same directory as the spectra.

20230216

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt

# local module, must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


# Create file list. Spectra in a (sub)folder.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")

# Line selection
linie, wellenlaenge = Linienlisten.linienauswahl()

# JD range of the series
flux, header = fits.getdata(filelist[0], header=True)
print('JD of the first spectrum', header['JD'])
flux, header = fits.getdata(filelist[-1], header=True)
print('JD of the last spectrum', header['JD'])

# Enter the range of observation times as JD
JD_anfang = float(input('Enter the start of the range of \
observation points (JD): '))
JD_ende = float(input('Enter the end of the range of \
observation times (JD): '))

bereich = float(input('Enter the speed range to be displayed around \
the selected line in km/s: '))

# Enter the system speed:
systemgeschwindigkeit = float(
    input("Enter the system speed in km/s: "))
systemgeschwindigkeit = systemgeschwindigkeit / 299772 * wellenlaenge

# Parameter for the distance between the spectra in the offset plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")


k = -1
# Work through the spectra list:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']
    JD_float = float(JD)

    if JD_float >= JD_anfang and JD_float <= JD_ende:
        k += 1
        wave = np.zeros(header['NAXIS1'])
        for j in range(header['NAXIS1']):
            wave[j] = lambda0 + j * step

        index_wellenlaenge = int(
            (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

        vstep = step / wellenlaenge * 299772
        pixbereich = int(bereich / vstep)

        fluxbereich = flux[(index_wellenlaenge - pixbereich)                           :(index_wellenlaenge + pixbereich)]
        wellenlaengenbereich = wave[(
            index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
        geschwindigkeitsbereich = (wellenlaengenbereich - wellenlaenge) / \
            wellenlaenge * 299772

        plt.figure(1, figsize=(7, 10))  # Overplot
        plt.plot(geschwindigkeitsbereich, fluxbereich, linewidth=.5)
        plt.xlabel('Velocity relative to the rest wavelength in km/s')
        plt.ylabel('relative intensity')
        plt.title(obj + ', ' + linie)
        plt.pause(.01)

        plt.figure(2, figsize=(7, 10))  # Overplot mit offset
        plt.plot(geschwindigkeitsbereich, fluxbereich +
                 k * offset, "-", label=JD, linewidth=1)
        # Beschriftung der einzelnen Spektren:
        plt.text(geschwindigkeitsbereich[1], 1.0 +
                 k * offset, format(JD, '.2f'), ha="left", size=7)
        plt.xlabel('Velocity relative to the rest wavelength in km/s')
        plt.ylabel('relative intensity')
        plt.title(obj + ', ' + linie)
        plt.axvline(x=0, color='k', linewidth=.05)
        plt.grid(visible=True, axis='x')
        plt.pause(.01)


frage = input(
    'Would you like to save the graphics as pdf and png? If yes, enter "y": ')
if frage == 'y':
    # plt.figure(1)
    # plt.savefig('Overplot_' + obj + '_' + linie + '_JD_range_' + str(JD_anfang) +
    #             '_' + str(JD_ende) + '.pdf')
    # plt.savefig('Overplot_' + obj + '_' + linie + '_JD_range_' + str(JD_anfang) +
    #             '_' + str(JD_ende) + '.png')
    plt.figure(2)
    # plt.savefig('OverplotWithOffset' + obj + '_'+linie + '_JD_range_' +
    #             str(JD_anfang) + '_' + str(JD_ende) + '.pdf')
    plt.savefig('OverplotWithOffset' + obj + '_'+linie + '_JD_range_' +
                str(JD_anfang) + '_' + str(JD_ende) + '.png')

plt.close('all')

print('\nEnd of program')
