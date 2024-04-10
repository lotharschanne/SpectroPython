#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates the equivalent distance of a line divided into areas for a time
series normalized to the continuum in fits format.
Input of the integration limits (in Angstroms) via a list which is to be
adapted to the respective case in line 24.
The EW calculation assumes that the same wavelength interval can be used to
calculate the integral for all spectra in the series, i.e. that no significant
RV changes occur. It is therefore best to use barycentrically corrected spectra.
The calculated EWs are saved in an asci file (csv format).

Created 2023-04-26

@author: lothar
"""

import numpy as np
from astropy.io import fits
import glob

###############################################################################
# List of wavelengths that serve as range limits for EW integration
WL = [6500.1, 6501.2, 6502.3]  # please adjust
##############################################################################

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
# print("\nList of spectra: \n")
# print(filelist)
print("Number of spectra: ", len(filelist), "\n")

fileobj = open('EWs.csv', 'w')
fileobj.write('List of range boundaries in Angstrom:,')
fileobj.write(str(WL)+'\n')
fileobj.write('Name of spectrum,')
fileobj.write('JD,')
fileobj.write('Equivalentwidth in Angström\n')


# Processing the filelist
for k in np.arange(len(filelist)):
    file = filelist[k]
    sp = fits.open(filelist[k], ignore_missing_end=True)
    try:
        JD = sp[0].header['JD']
    except:
        JD = None

    fileobj.write(file)
    fileobj.write(',')
    fileobj.write(str(JD))
    fileobj.write(',')
    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )

    # Close the fits-file:
    sp.close()

    # Processing the integration ranges:
    for m in np.arange(len(WL)-1):
        EW = 0
        for n in np.arange(len(wave)):
            if wave[n] >= WL[m] and wave[n] <= WL[m+1]:
                EW += (1 - flux[n])*sp[0].header['CDELT1']
        fileobj.write(str(EW))
        fileobj.write(',')
    fileobj.write('\n')

fileobj.close()
