#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d spectra in fits format.
They must all have the same
header['CRVAL1'], header['CDELT1'], header['NAXIS1'], header['CRPIX1'].

These variables are printed for each spectrum (control). With the
input of 'y', all fluxes of the spectrum series are averaged and the mean
spectrum is saved as a fit.

Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input('Path and name of the spectra (use wildcards): ')
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
print('\List of spectra: \n')
print(filelist)
print('\nNumber of spectra: ', len(filelist), '\n')


print('Header entries')
print('CRVAL1', 'CDELT1', 'NAXIS1', 'CRPIX1')

for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)

    print(header['CRVAL1'], header['CDELT1'],
          header['NAXIS1'], header['CRPIX1'])

frage = input('Should the spectra be averaged? Then enter y: ')

if frage == 'y':
    fluxsumme = np.zeros(int(header['NAXIS1']))
    for i in range(len(filelist)):
        flux, header = fits.getdata(
            filelist[i], header=True, ignore_missing_end=True)
        fluxsumme += flux
    fluxsumme = fluxsumme / len(filelist)
    fits.writeto('mean_'+files, fluxsumme, header, overwrite=True,output_verify='silentfix')
