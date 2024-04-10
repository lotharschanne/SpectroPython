#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in a 1d spectrum in the form of a fits file.

The spectrum can be convolved with a standard deviation, expressed by the
resolution R. The convolved fits file is saved.

Stand 20220105
@author: Lothar Schanne
"""
# from __future__ import print_function, division
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii, fits

spectrum_name = input('Enter path and file name: ')

R = float(input('Desired R: '))

# Read spectrum
sp = fits.open(spectrum_name)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')
std = sp[0].header['CRVAL1']/R/sp[0].header['CDELT1']

#   Convolve with astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
sp[0].data = convolve(sp[0].data, kernel, normalize_kernel=True,
                      boundary='extend')

name = spectrum_name.rstrip('.fits')[0]+'_convolved_R'+str(R)+'.fits'


# Speichern des convolved
fits.writeto(name, sp[0].data, sp[0].header,
             overwrite=True, output_verify='silentfix')
