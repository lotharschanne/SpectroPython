#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in a series of heliocentrically corrected spectra in fits format,
Input of an evaluation period as JD start and JD end.
Display of several selectable lines in velocity space in one plot,
plotted on top of each other with a selectable offset.
The plots can be saved as pdf and png.

The file ListOfSpectralLines.py must be in the same directory as the spectra.

20230221

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


#  Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print("Number of spectra: ", len(filelist), "\n")


# JD range of the series
flux, header = fits.getdata(filelist[0], header=True)
print('JD of the first spectrum', header['JD'])
flux, header = fits.getdata(filelist[-1], header=True)
print('JD of the last spectrum', header['JD'])

# Enter the range of observation times as JD
JD_anfang = float(input('Enter the start of the range of \
observation times (JD): '))
JD_ende = float(input('Enter the end of the range of \
observation times (JD): '))

wahl = True
linie = []
wellenlaengen = []
while wahl:
    # Selecting a line
    linie, wellenlaenge = Linienlisten.linienauswahl()
    wellenlaengen.append(wellenlaenge)
    wahl = bool(input('Select another line? Then enter y,\
otherwise press return: '))

Farbenwahl = ['k', 'r', 'b', 'g', 'c', 'm', 'y']
Farben = Farbenwahl[0:len(wellenlaengen)]

bereich = float(input('Enter the speed range to be displayed around \
the selected line in km/s: '))

# Enter the system speed:
systemgeschwindigkeit = float(
    input("Enter the system speed in km/s: "))

# Parameter for the distance between the spectra in the offset plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

k = -1  # Initialization of the counter for the spectra in the JD range
# Process the spectra list:
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

        for wert in wellenlaengen:
            index_wellenlaenge = int(
                (wert - lambda0 +
                 (systemgeschwindigkeit / 299772 * wert)) / step)

            vstep = step / wert * 299772
            pixbereich = int(bereich / vstep)

            fluxbereich = flux[(index_wellenlaenge - pixbereich):
                               (index_wellenlaenge + pixbereich)]
            wellenlaengenbereich = wave[(
                index_wellenlaenge - pixbereich):(index_wellenlaenge +
                                                  pixbereich)]
            geschwindigkeitsbereich = (wellenlaengenbereich - wert)\
                / wert * 299772

            plt.figure(1, figsize=(7, 10))  # Overplot with offset
            if k == 0:
                plt.plot(geschwindigkeitsbereich, fluxbereich +
                         k * offset, "-", c=Farben[wellenlaengen.index(wert)],
                         linewidth=1, label=str(wert))
            else:
                plt.plot(geschwindigkeitsbereich, fluxbereich +
                         k * offset, "-", c=Farben[wellenlaengen.index(wert)],
                         linewidth=1)
            # Labeling of the individual spectra:
            if wert == wellenlaengen[-1]:
                plt.text(geschwindigkeitsbereich[1], 1.0 +
                         k * offset, format(JD, '.2f'), ha="left", size=7)
                plt.xlabel(
                    'Velocity relative to the rest wavelength in km/s')
                plt.ylabel('relative intensity')
                plt.title(obj + ', ' + str(wellenlaenge))
                plt.axvline(x=0, color='k', linewidth=.05)
                plt.grid(visible=True, axis='x')
            plt.pause(.01)
        plt.legend()


frage = input(
    'Would you like to save the graphic as pdf and png? If yes, enter "y": ')
if frage == 'y':
    plt.savefig('OverplotWithOffset_' + obj + '_' + str(wellenlaenge) +
                '_JD_range_' + str(JD_anfang) + '_' + str(JD_ende) + '.pdf')
    plt.savefig('OverplotWithOffset_' + obj + '_'+str(wellenlaenge) +
                '_JD_range_' + str(JD_anfang) + '_' + str(JD_ende) + '.png')

plt.close('all')

print('\nEnd of program')
