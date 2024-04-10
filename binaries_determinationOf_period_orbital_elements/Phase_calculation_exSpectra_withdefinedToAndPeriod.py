#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script reads in a catalog of fits-1d spectra.
Using the parameters specified in the script (which may need to be adjusted)
(period and T0) are then used to convolve the observation times, and the results
(phase, JD) are saved in an ascii file (csv format).
The spectra must contain a header entry called 'JD'.

Created on 20231219

@author: lothar
"""
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import numpy as np
import glob

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

Periode = 680.168  # adjust   !!!!!!!!!!!!!!!!!!!!!!
T0 = 2454433.2     # adjust   !!!!!!!!!!!!!!!!!!!!!!

# Create file list. All spectra in one (sub)folder.
files = input('Path and name of the fits-files (use wildcards) : ')
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print('\nList of spectra: \n')
print('Number of spectra: ', len(filelist), '\n')
jd = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    try:
        jd[i] = header['JD']
    except:
        print('Spectrum '+filelist[i]+'has no JD entry')
        jd[i] = None
        pass

# Folding (Phase)
phases = pyasl.foldAt(jd, Periode, T0)


data = pd.DataFrame({'Spectrum': filelist, 'JD': jd, 'Phase': phases})
datafile = data.to_csv('Phase_JD_exSpectra_perToAndPeriod.csv')
