#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates for a series of 1d-fits spectra using the beta*sigma method
the SNR of each spectrum. Saves the SNR in a file (comma separated, .csv).


Stand 20221105

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits, ascii
import glob

# Create file list. Spectra in a (sub)folder.
files = input('Path and name of the spectra (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print('\nSpectrum list: \n')
print('Number of spectra: ', len(filelist), '\n')

fluxmean = np.zeros(len(filelist))
nstd = np.zeros(len(filelist))
nstdstd = np.zeros(len(filelist))
snr = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    # Estimate noise using robust estimate
    beq = pyasl.BSEqSamp()
    # Define order of approximation (use larger values such as 2,3, or 4 for
    # faster varying or less well sampled data sets; also 0 is a valid order)
    N = 1
    # Define 'jump parameter' (use larger values such as 2,3, or 4 if correlation
    # between adjacent data point is suspected)
    j = 3
    # Estimate noise assuming equidistant sampling (often a good approximation even
    # if data are not strictly equidistant) and robust estimation (often advantageous
    # in working with real data)
    nstd[i], nstdstd[i] = beq.betaSigma(flux, N, j, returnMAD=True)

    fluxmean[i] = flux.mean()
    snr[i] = fluxmean[i] / nstd[i]

    print('\nSpectrum: ', filelist[i])
    print("Estimated noise std = %7.5f +/- %7.5f" % (nstd[i], nstdstd[i]))
    print('SNR = ', int(snr[i]))
    print('used Parameters N = ', N, ' and j = ', j)

ascii.write([filelist, fluxmean, nstd, nstdstd, snr], files+'_SNR'+'.csv', overwrite=True,
            names=['Spectrum', 'mean Flux', 'noise', 'noise_std', 'SNR'], format='csv')
