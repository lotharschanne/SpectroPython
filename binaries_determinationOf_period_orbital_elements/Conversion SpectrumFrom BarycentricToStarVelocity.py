#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Required:
    Spectra as fits
    Table RVs and phases of the star calculated with the orbit parameters of
    the orbit

@author: lothar
"""


import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the fits difference spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
# print("\nList of spectra: \n")
# print(filelist, end='\n')
print("Number of Spectra: ", len(filelist), "\n")

tabelle = ascii.read(
    'theoretischeRVsAlgolABCzuJDZeitpunktenAusSpektren.csv', format='csv')  # adjust !!!!
print('Column names of the table:')
print(tabelle.colnames)


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    for j in range(len(tabelle)):
        if tabelle['JD'][j] == header['JD']:
            # if another RV is desired adjust !!!
            rv_A = tabelle['RV_A'][j]
            wave = np.zeros(header["NAXIS1"], dtype=float)
            for k in range(len(wave)):
                if 'CRPIX1' not in header:
                    header['CRPIX1'] = 1
                wave[k] = header["CRVAL1"] + \
                    (k + 1 - header["CRPIX1"]) * header["CDELT1"]
            # Shift the spectrum
            flux_rv, wave_rv = pyasl.dopplerShift(
                wave, flux, -rv_A, edgeHandling="firstlast")
            header["CRVAL1"] = wave[0]
            header['CRPIX1'] = 1
            fits.writeto(filelist[i].rsplit('.fit')[0]+'on_star_related.fit',
                         flux_rv, header, overwrite=True,
                         output_verify="silentfix")  # Customize file name !!!
            ascii.write([wave, flux_rv],
                        filelist[i].rsplit('.fit')[0]+'on_star_related.csv',
                        overwrite=True,
                        names=['WAVE', 'FLUX'],
                        format="csv")  # Customize file name  !!!
        else:
            pass
