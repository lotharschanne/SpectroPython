#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import of a 1d spectrum in the form of an ascii table with the extension .dat,
2-column with the column headings 'WAVE' and 'FLUX'.
The spectrum can be convolved with a standard deviation
(stddev in the unit of the step).
The convolved flux is stored together with the wavelength in a two-column
ascii table (spacer = tab). With the column headings
'WAVE' and 'FLUX'.

Stand 20180824
@author: Lothar Schanne
"""
# from __future__ import print_function, division
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii

spectrum_name = input('Enter path and file name: ')

std = float(input('Standard deviation for the convolution in multiples of\
 increment: '))

# Read spectrum
spectrum = ascii.read(spectrum_name, guess=True)

#   Convolve with astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
convoluted = convolve(spectrum['FLUX'], kernel, normalize_kernel=True,
                      boundary='extend')

name = spectrum_name.rstrip('.dat')[0]

# Grafik
fig = plt.figure(1, figsize=(14, 10))
plt.suptitle('Spectrum '+spectrum_name)
plt.subplot(2, 1, 1)
plt.plot(spectrum['WAVE'], spectrum['FLUX'])
plt.xlabel('Wavelength in angstroms')
plt.ylabel('relative Flux')
plt.subplot(2, 1, 2)
plt.plot(spectrum['WAVE'], convoluted)
plt.xlabel('Wavelength in angstroms')
plt.ylabel('relative Flux')
plt.title(name)
# saving
plt.savefig(name+'.png')
# plt.savefig(name+'.pdf')
plt.pause(.1)

# Save the convolved
ascii.write([spectrum['WAVE'], convoluted],
            name+'_convolved.dat',
            overwrite=True,
            names=['WAVE', 'FLUX'],
            format='tab')
