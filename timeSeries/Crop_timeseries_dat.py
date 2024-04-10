#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create file list for wavelength calibrated 1d spectra of a time series in the
ascii format (.dat). Plot the spectra and determine the common wavelength range.
Trimming of all spectra to a selectable wavelength range and
saving all clipped spectra as fits with wavelength range in the file name

20220125
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

# Create file list.
files = input("Path and name of the dat spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
# print("\nList of spectra: \n")
# print(filelist)
print("Number of spectra: ", len(filelist), "\n")
print('Please wait. Calculations are running.')


fig = plt.figure(figsize=(14, 20))
# Einlesen von flux und header:

tab = ascii.read(filelist[0], format='tab')
print('Table column names: ', tab.colnames)
flux = tab['FLUX']
wave = tab['WAVE']
plt.plot(wave, flux, linewidth=1)
print(filelist[0], wave[0], wave[-1])
print('Please wait. Calculations are running.')

plt.grid(True)
plt.xlabel("Wavelength [Angstrom]")

plt.pause(20)  # During this time, the graphic can be actively edited.

# Generate the spectrum section, plot and save as fit
min = 0
max = 10000
for i in range(len(filelist)):
    wave = ascii.read(filelist[i])['WAVE']
    minimum = wave[0]
    maximum = wave[-1]
    if minimum > min:
        pass
    else:
        min = minimum
    if maximum < max:
        pass
    else:
        max = maximum

print("\nCommon wavelength range: ",
      int(min), ' to ', int(max),  ' Angstrom')

print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Begin: "))
b = float(input("End: "))

plt.close()

for i in range(len(filelist)):
    tab = ascii.read(filelist[i], format='tab')
    print('Table column names: ', tab.colnames)
    flux = tab['FLUX']
    wave = tab['WAVE']
    for j in range(len(wave)):
        if wave[j] >= a:
            aindex = j
            break
    for j in range(len(wave)):
        if wave[j] >= b:
            bindex = j
            break

    newflux = flux[aindex:bindex]
    newwave = wave[aindex:bindex]

    ascii.write([newwave, newflux],
                filelist[i].rstrip('.dat')+'_'+str(a)+'_'+str(b)+'.dat',
                overwrite=True,
                names=["WAVE", "FLUX"],
                format="tab",)

    plt.figure()
    plt.plot(newwave, newflux, 'k-', linewidth=1)
    plt.title(filelist[i], fontsize=9)
    plt.xlabel('Wavelength [Angstrom]')
    plt.ylabel('relative Flux')
    plt.ylim(-.2, .2)  # adjust !!!!!!!!!!!!!!!!!!
    plt.grid(True)
    plt.savefig(filelist[i].rstrip('.dat')+'_'+str(a)+'_'+str(b)+'.png',
                format='png')
    plt.pause(.5)
    plt.close()


print('End of the program')
plt.close('all')
