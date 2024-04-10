#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Import of a series of 1d-fits files and automatic normalization to the mean value of the
average value of the fluxes of the spectrum. Saving the calculated fits with the
name suffix _1Pnorm.

20240408

@author: lothar
"""

from astropy.io import fits
import glob


# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the files (use wildcards) : ')
filelist = glob.glob(files)

print('Number of spectra: ', len(filelist), '\n')


# Read header and flux
for i in range(len(filelist)):
    print(filelist[i], ':')
    flux, header = fits.getdata(filelist[i], header=True)
    flux = flux / flux.mean()
    newfile = filelist[i].rstrip('.fits') + '_1Pnorm' +  '.fits'
    fits.writeto(newfile, flux, header,
                 overwrite=True, output_verify='silentfix')
