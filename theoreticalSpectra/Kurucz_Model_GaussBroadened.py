#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Loading a Kuruzc star atmosphere model (1d spectrum) after selecting Teff and
logg (internet connection required).
Convolution with a Gaussian function to spectrograph resolution,
selection of the wavelength range and adjustment to a desired step width
(angstrom/pixel).
Saving as ascii-tab-table and as .fit with the filename to be entered.
Saving the plot as pdf.

Created on Fri Apr  2 13:53:54 2021

@author: lothar
"""
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from astropy.io import ascii, fits


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

# Auswahl des Modells mittels Teff, logg und Metallizität
model = pyasl.resBased.spectralLib.SpectralLib(lib='A')
print('List of available models (parameters and file names)')
model.listInventory()

Teff = float(input('Enter Teff: '))
logg = float(input('Enter logg: '))
met = 0.0

Model = model.requestModel(Teff, logg, met)
wave, flux = model.read1dFitsSpec(Model)

specname = input('Enter file name for saving: ')+'_'

print('Step width of the model spectrum = ', wave[1]-wave[0])
R = float(input('Enter the desired relative resolution R (e.g. 10000): '))
newstep = float(
    input('Enter the desired step width of the final result: '))
spectrum_name = specname + str(Teff) + '_' + str(logg) + '_R' + str(int(R))

print('\nEnter the desired wavelength range in angstroms:')
a = float(input('Begin: '))
b = float(input('End: '))

for m in range(len(wave)):
    if wave[m] > a:
        aind = m
        break
for m in range(len(wave)):
    if wave[m] > b:
        bind = m
        break
wave = wave[aind:bind]
flux = flux[aind:bind]


convolutedFlux, fwhm = pyasl.instrBroadGaussFast(
    wave, flux, R, edgeHandling='firstlast', fullout=True, maxsig=5.0)
print('FWHM of the Gaussian broadened spectrum = ', fwhm)

data, dt = pyasl.binningx0dt(
    wave, convolutedFlux, x0=wave[0], dt=newstep, useBinCenter=True, useMeanX=True)
name = spectrum_name + '_convolved.dat'
ascii.write([data[:, 0], data[:, 1]], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
hdu = fits.PrimaryHDU()
hdu.header['CRVAL1'] = data[0, 0]
hdu.header['CRPIX1'] = 1
hdu.header['NAXIS1'] = len(data)
hdu.header['CDELT1'] = newstep
hdu.header['NAXIS'] = 1
name = spectrum_name + '_convolved.fit'
fits.writeto(name, data[:, 1], hdu.header,
             overwrite=True, output_verify='silentfix')


# Grafik
fig, ax = plt.subplots(nrows=2, ncols=1)
fig.suptitle('Spectrum: ' + spectrum_name, fontsize=12)

ax[0].plot(wave, flux, linewidth=.3)
ax[0].set(ylabel='relative Flux')
ax[0].set_title('Kurucz-Model', fontsize=8)

ax[1].plot(data[:, 0], data[:, 1], linewidth=.3)
ax[1].set(xlabel='Wavelength in angstroms', ylabel='relative Flux')
ax[1].set_title('convolved', fontsize=8)
# ax[1].set_xlim(6500,6700)

plt.savefig(spectrum_name + '_convolved.pdf')
