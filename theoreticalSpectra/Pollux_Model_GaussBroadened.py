#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loading a Pollux model spectrum (format .spec).
Convolution with a Gaussian function to spectrograph resolution,
selection of the wavelength range and
adjustment to a desired step width (Angstrom/pixel).
Saving as ascii-tab-table and as .fit with the filename to be entered.
Saving the plot as pdf.

Created on Fri Apr  2 13:53:54 2021

@author: lothar
"""
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
import numpy as np


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

file = input("Enter the path and file name of the Pollux spectrum: ")
table = ascii.read(file)

specname = input('Enter the file name for saving the generated spectrum: ')+'_'

print('\nEnter the desired wavelength range in angstroms:')
a = float(input('Begin: '))
b = float(input('End: '))

newtable = table[table['col1'] >= a]
newtable = newtable[newtable['col1'] <= b]
wave = newtable['col1']
flux = newtable['col3']

print('Step size of the model spectrum = ', wave[1]-wave[0])

R = float(input('Enter the desired relative resolution R (e.g. 10000): '))
newstep = float(
    input('Enter the desired step width in angstroms of the final result: '))
spectrum_name = specname + '_R' + str(int(R))

convolutedFlux, fwhm = pyasl.instrBroadGaussFast(
    wave, flux, R, edgeHandling='firstlast', fullout=True, maxsig=5.0, equid=True)
print('FWHM of the Gaussian broadened spectrum = ', fwhm)

neue_wellenlaengen = np.arange(wave[0], wave[-1], newstep)
newwave, newflux = pyasl.equidistantInterpolation(wave, convolutedFlux, neue_wellenlaengen)

name = spectrum_name + '_convolved.dat'
ascii.write([newwave, newflux], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
hdu = fits.PrimaryHDU()
hdu.header['CRVAL1'] = newwave[0]
hdu.header['CRPIX1'] = 1
hdu.header['NAXIS1'] = len(newwave)
hdu.header['CDELT1'] = newstep
hdu.header['NAXIS'] = 1
name = spectrum_name + '_convolved.fit'
fits.writeto(name, newflux, hdu.header,
              overwrite=True, output_verify='silentfix')


# Grafik
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=True)
fig.suptitle('Spectrum: ' + file, fontsize=10)

ax[0].plot(wave, flux, linewidth=.5)
ax[0].set(ylabel='relativer Flux')
ax[0].set_title('Pollux-Model', fontsize=8)

ax[1].plot(newwave, newflux, linewidth=.5)
ax[1].set(xlabel='Wavelength in angstroms', ylabel='relative Flux')
ax[1].set_title('convolved to R ' + str(R), fontsize=8)
# ax[1].set_xlim(6500,6700)
plt.pause(1)

plt.savefig(spectrum_name + '_convolved.png', format='png')
# plt.savefig(spectrum_name + '_convolved.pdf', format='pdf')
