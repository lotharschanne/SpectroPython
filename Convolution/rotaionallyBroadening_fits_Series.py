#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The spectra are rotationally broadened by entering the limb-darkening
coefficient and the rotation speed.
The convolved fits files are saved, with the term 'rotationally_broadened'
added to the file names.

Stand 20220105
@author: Lothar Schanne
"""

import matplotlib.pyplot as plt
from astropy.io import fits
import glob
import numpy as np
from PyAstronomy import pyasl

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

epsilon = float(
    input('Enter the limb-darkening coefficient, between 0 and 1: '))
vsini = float(
    input('Enter the projected rotational speed vsini in km/s: '))

# Ausdruck der Liste
print("\nSpectra list: \n")
print("Number of spectra: ", len(filelist), "\n")

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    step = header['CDELT1']
    if 'CRPIX1' in header:
        refpix = header['CRPIX1']
    else:
        header['CRPIX1'] = 1
        refpix = 1
    wave_erstesPix = header['CRVAL1'] - step*(refpix - 1)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = wave_erstesPix + k * step

    # Rotational broadening
    flux_broadened = pyasl.rotBroad(wave, flux, epsilon, vsini)

    name = filelist[i].rstrip('.fits')[0] + '_rotationally_broadened' + '.fits'

    # Saving the convolved
    fits.writeto(name, flux_broadened, header,
                 overwrite=True, output_verify='silentfix')
