#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert a series of wavelength calibrated 1d spectra in fits format into text
format .csv (comma-separated) and .dat format (tab-separated).
With column headers 'WAVE' and 'FLUX'.

Stand 20180815
@author: lothar
"""

import numpy as np
from astropy.io import fits
from astropy.io import ascii
import glob

# Create file list. Spectra in a (sub)folder.
files = input('Path and name of the spectra (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Ausdruck der Liste
print('\List of spectra: \n')
print(filelist)
print('Number of spectra: ', len(filelist), '\n')
print('Please wait. Calculations are running.')

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    #   Check for the necessary header entries
    print('\nSpektrum: ', filelist[i])
    print('Output of the header entries required for the wavelength calculation:')
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
        print('AStarting wavelength, CRVAL1:             ', crval)
    else:
        print('CRVAL1 is missing in the header !')
    if 'CDELT1' in header:
        cdel = header['CDELT1']
        print('Step width of the wavelength, CDELT1:   ', cdel)
    else:
        print('CDELT1 is missing in the header !')
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1

    #   Creation of numpy arrays with the wavelengths and fluxes of the spectrum
    wave = np.ones(nax, dtype=float)
    for j in range(nax):
        wave[j] = crval + (j - header['CRPIX1'] + 1) * cdel
    ascii.write([wave, flux], filelist[i].strip('.fit')+'.csv', overwrite=True,
                names=['WAVE', 'FLUX'], format='csv')
    ascii.write([wave, flux], filelist[i].strip('.fit')+'.dat', overwrite=True,
                names=['WAVE', 'FLUX'], format='tab')

print('End of the program')
