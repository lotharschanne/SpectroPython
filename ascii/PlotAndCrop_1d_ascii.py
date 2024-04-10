#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PlotAndCrop_1d_ascii.py

The script reads an ascii file and plots it. The ascii file has 2 columns
with the column headings 'WAVE' and 'FLUX'. The extension of the file is .dat.
The separator is determined automatically with guess=True.

Optionally, a section (wavelength range) can be created.
This is saved as .dat and plotted.

Created on Sun Jul 29 15:43:51 2018
@author: lothar
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


# Input of the ascii file
# This must contain 2 arbitrarily but uniformly separated columns WAVE and FLUX
spectrum_name = input('Enter the path and name of the ascii file: ')
spectrum = ascii.read(spectrum_name, guess=True)

fig = plt.figure(figsize=(14, 10))
plt.subplot(2, 1, 1)
plt.plot(spectrum['WAVE'], spectrum['FLUX'])
plt.xlabel('Wavelength [Angstroem]')
plt.ylabel('ADU')
plt.title(spectrum_name)
plt.grid(True)
plt.show(block=True)


# Create the section, plot and save as ascii-file
antwort = input('Would you like to create a section? y/n:  ')
if antwort == 'y':
    print('Specification of the wavelength range to be transferred: ')
    aindex = int(input('Begin of the spectrum: '))
    bindex = int(input('End of the spectrum: '))
    filename = spectrum_name.rstrip('.dat')+'_'+str(aindex)+'_'+str(bindex)
    # Search indices
    file_wave = np.ones(len(spectrum))
    file_flux = np.ones(len(spectrum))
    for m in range(len(spectrum)):
        file_wave[m] = spectrum['WAVE'][m]
    for m in range(len(spectrum)):
        file_flux[m] = spectrum['FLUX'][m]
    for m in range(len(spectrum)):
        if file_wave[m] > aindex:
            aind = m
            break
    for m in range(len(spectrum)):
        if file_wave[m] > bindex:
            bind = m
            break
    file_wave = file_wave[aind:bind]
    file_flux = file_flux[aind:bind]
    ascii.write([file_wave, file_flux], filename+'.dat', overwrite=True,
                names=['WAVE', 'FLUX'], format='tab')
    plt.subplot(2, 1, 2)
    plt.plot(file_wave, file_flux)
    plt.xlabel('Wavelength [Angstroem]')
    plt.ylabel('ADU')
    plt.title(filename)
    plt.grid(True)

    # fig.savefig(filename+'.pdf')
    fig.savefig(filename+'.png')

    plt.pause(.1)
