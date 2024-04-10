#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards): ")
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
print("\nNumber of spectra: ", len(filelist), "\n")
print("\nSpectrum list: \n")
print(filelist)


date_obs = list(range(len(filelist)))
jd = list(range(len(filelist)))

print("\nHeader entries")

# Read in and output observation time
print("\nSpectrum", "DATE-OBS", 'JD')

for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)
    if "DATE-OBS" in header:
        date_obs[i] = header["DATE-OBS"]
    elif "DATE_OBS" in header:
        date_obs[i] = header["DATE_OBS"]
    elif "MJD" in header:
        date_obs[i] = header["MJD"]
    elif "JD-MID" in header:
        date_obs[i] = header["JD-MID"]
    else:
        print("No observation time in the header")

    if "JD" in header:
        jd[i] = header['JD']
    else:
        print("No JD in the header")

    print(filelist[i], date_obs[i], jd[i])

# Calculate and output the spectral range
print("    Spectrum file    ", "Begin   ", "End")
for i in range(len(filelist)):
    flux, header = fits.getdata(
        filelist[i], header=True, ignore_missing_end=True)

    Anfang = header["CRVAL1"]
    Ende = header["CRVAL1"] + header["NAXIS1"] * header["CDELT1"]
    print(filelist[i], f"{Anfang:.1f}, {Ende:.1f}")

# Saving as an ascii file
ascii.write(
    [filelist, date_obs, jd],
    "ObservationDatesAndJD" + ".dat",
    overwrite=True,
    names=["Spectrum", "DATE-OBS", 'JD'],
    format="tab",
)
