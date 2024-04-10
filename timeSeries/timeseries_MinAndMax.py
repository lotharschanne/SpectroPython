#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Reads in a time series of 1d spectra in fits format.
Calculates the wavelength and the flux of the minimum and maximum in the
spectrum.

Stand 20231214
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits, ascii
import glob


# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the file (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

minimumwave = np.zeros(len(filelist))
maximumwave = np.zeros(len(filelist))
minimumflux = np.zeros(len(filelist))
maximumflux = np.zeros(len(filelist))
JD = np.zeros(len(filelist))

# Read header and flux
for i in range(len(filelist)):
    print(i, filelist[i])
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections and save them as fit:
    step = header["CDELT1"]
    refpix = header["CRPIX1"]
    # JD of the observation, used to label the spectra:
    if 'JD' in header:
        JD[i] = header['JD']
    elif 'MJD-OBS' in header:
        JD[i] = header['MJD-OBS']
    elif 'DATE-OBS' in header:
        JD[i] = header['DATE-OBS']
    else:
        print('\nNo observation date in header')

    wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)

    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = wave_erstesPix + k * step

    minimumflux[i] = flux.min()
    maximumflux[i] = flux.max()
    minimumwave[i] = wave[flux.argmin()]
    maximumwave[i] = wave[flux.argmax()]


# Save as an ascii file
ascii.write(
    [filelist, JD, minimumwave, minimumflux, maximumwave, maximumflux],
    "MinMax_atJD" + ".dat",
    overwrite=True,
    names=["Spectrum", "JD", 'Wavelength minimum', 'Flux minimum',
           'Wavelength maximum', 'Flux Maximum'],
    format="tab",
)
