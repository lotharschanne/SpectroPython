#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wellenlaengenbereich_Zeitserie_fits.py

Create file list for wavelength calibrated 1d spectra of a time series in the
fits format. Save the filelist with information on the end and start of the
wavelength scale.
Print out the common wavelength range.
20220125
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) :")
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print(filelist)
print("Number of spectra: ", len(filelist), "\n")

begin = np.zeros(len(filelist))
end = np.zeros(len(filelist))
print('Spectrum          begin   end')

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    begin[i] = header["CRVAL1"]
    end[i] = header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"]
    print(filelist[i], begin[i].round(2), end[i].round(2))

ascii.write(
    [filelist, begin.round(2), end.round(2)],
    "Wavelength_ranges.dat",
    names=["Spectrum", "Begin", "End"],
    format="tab",
    overwrite=True,
)

min = 0
max = 10000
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if header["CRVAL1"] > min:
        min = header["CRVAL1"]
    if header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"] < max:
        max = header["CRVAL1"] + header["CDELT1"] * header["NAXIS1"]
print("\nGCommon wavelength range: ", min, max)
