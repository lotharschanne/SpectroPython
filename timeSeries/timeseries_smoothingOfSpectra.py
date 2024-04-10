#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read in a series of 1d-fits files and smooth the flux with a certain window
width. Saving the smoothed fits.

Created on Fri Dec 10 16:10:46 2021

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the file (use wildcards) : ')
filelist = glob.glob(files)

obj = input('Enter the name of the object:')
namenszusatz = input('Enter the name suffix for the spectra to be saved: ')

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')

# The smoothing window must be an odd number, adjust !!!!!!!!!!!!!!!!!!!
Fenster = 21

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i], ':')
    flux, header = fits.getdata(filelist[i], header=True)
    flux_smoothed = pyasl.smooth(flux, Fenster, 'flat')
    newfile = obj + '_' + header['JD'] + '_' + namenszusatz + '.fits'
    fits.writeto(newfile, flux_smoothed, header,
                 overwrite=True, output_verify='silentfix')
