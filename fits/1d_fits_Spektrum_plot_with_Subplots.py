#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Viewing a wavelength-calibrated 1d spectrum, reading out and displaying of
the header data.
Plotting of the entire spectrum and a division into a selectable number of
sections.

Created on Sunday Jul 22 2018

@author: lothar
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Create filelist
# Path and name of the spectra, please customize files
file = input('Enter path and file name: ')

#   Importing headers and data
flux, header = fits.getdata(file, header=True)

print('Minimum and maximum in flux: ', flux.min(), '  ', flux.max())

#   Header-Check Spectrum
header_flag = input('Would you like to see the header in full? y/n:')
if header_flag == 'y':
    print('Header data: ')
    print(header)

#   Check for the necessary header entries
print('\nOutput of the header entries required for the wavelength calculation:')
if 'NAXIS' in header:
    print('Dimension, NAXIS:                        ', header['NAXIS'])
else:
    print('This is not a 1d spectrum!')

if 'NAXIS1' in header:
    nax = header['NAXIS1']
    print('Number of values (abscissa), NAXIS1:     ', nax)
else:
    print('NAXIS1 is missing in the header !')

if 'CRVAL1' in header:
    crval = header['CRVAL1']
    print('Initial wavelength, CRVAL1:             ', crval)
else:
    print('CRVAL1 is missing in the header !')

if 'CDELT1' in header:
    cdel = header['CDELT1']
    print('SStep size of the wavelength, CDELT1:    ', cdel)
else:
    print('CDELT1 is missing in the header !')

if 'CRPIX1' not in header:
    header['CRPIX1'] = 1

#   Creating arrays with the wavelengths and fluxes of the spectrum
wave = np.ones(nax, dtype=float)
crval = crval + (1 - header['CRPIX1']) * cdel
for i in range(nax):
    wave[i] = crval + i * cdel
# The wave list contains the wavelengths of the pixels
# The corresponding intensities in the flux list


# Plot full spectrum
fig = plt.figure(1, figsize=(10, 10))
plt.plot(wave, flux)
plt.xlabel('Wavelength [Angström]')
plt.ylabel('ADU')
plt.title('Spectrum '+file)
plt.grid(True)
# fig.savefig(file.strip('.fit')+'.pdf', format='pdf')
fig.savefig(file.strip('.fit')+'.png', format='png')
plt.pause(1)

#   n Subplots
n = int(input('Enter the number of desired subplots: '))
if n > 1:
    pin = nax/n
    fig, ax = plt.subplots(n, figsize=(10, 2*n))
    fig.subplots_adjust(bottom=0.15, left=0.2)
    for i in range(n):
        links = int(i*pin)
        rechts = int((i+1)*pin)
        ax[i].plot(wave[links:rechts], flux[links:rechts])
        ax[i].grid(True, linestyle='-.')
        ax[i].set_ylim(0.95*flux.min(), 1.05*flux.max())
    ax[0].set_title('Spectrum '+file)
    ax[n-1].set_xlabel('Wavelength [Angström]')
    fig.savefig(file.strip('.fit')+'_zoom.pdf')
    plt.pause(1)

# Plot of a special wavelength range
frage = input('Would you like to zoom in on a specific wavelength range?\n\
The range can be selected based on the already saved pdf.y/n: ')
if frage == 'y':
    a = int(input('Beginning of the spectrum: '))
    b = int(input('End of the spectrum: '))
    aindex = int((a - crval) / cdel)
    bindex = int((b - crval) / cdel)
    fig = plt.figure(3, figsize=(10, 10))
    plt.plot(wave[aindex:bindex], flux[aindex:bindex])
    plt.xlabel('Wavelenght [Angström]')
    plt.ylabel('ADU')
    plt.title('Spectrum ')
    fig.savefig(file.rstrip('.fit')+'Cutout_1.png', format='png')

plt.pause(1)
