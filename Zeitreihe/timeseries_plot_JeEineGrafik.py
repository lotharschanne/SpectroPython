#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
timeseries_plot.py

Creates a file list for wavelength calibrated 1d_spectra (fits-format) of a time series.

Plots the spectra and save the graph.

Release 28.3.2021
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the file (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')


# Parameter für den Abstand zwischen den Spektren im zweiten Plot
# bitte anpassen
offset = 0.5

# Read header and flux
for i in range(len(filelist)):
    fig = plt.figure()
    print(filelist[i], ':')
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step = header['CDELT1']
    refpix = header['CRPIX1']
    wave_erstesPix = header['CRVAL1'] - step*(refpix - 1)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = wave_erstesPix + k * step

    plt.plot(wave, flux, '-', label=filelist[i], linewidth=.4)
    plt.grid(True)
    plt.title(filelist[i])
    plt.xlabel('Wavelength in Angström')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    # Legend can be adjusted, e.g. font size (fontsize)
    # plt.legend(bbox_to_anchor=(0., 1.02, 1., .402), loc=3,
    #            ncol=2, mode="expand", borderaxespad=1., fontsize=6)

    fig.savefig(filelist[i]+'.pdf', format = 'pdf')
