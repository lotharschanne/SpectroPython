#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create file list for wavelength calibrated 1d spectra of a time series in the
fits format. Plot the spectra and determine the common wavelength range.
Trimming of all spectra to a selectable wavelength range and saving all clipped
spectra as fits with wavelength range in the file name.

20220125
@author: lothar
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist)
print("Number of spectra: ", len(filelist), "\n")
print('Please wait. Calculations are running.')


fig = plt.figure(figsize=(14, 20))
# Read in flux and header:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    wave = np.zeros(header["NAXIS1"], dtype=float)
    for k in range(len(wave)):
        wave[k] = header["CRVAL1"] + \
            (k + 1 - header["CRPIX1"]) * header["CDELT1"]
    plt.plot(wave, flux, linewidth=1)
    print(filelist[i], wave[0])
print('Please wait. Calculations are running.')

# Customize plot properties (grid, labels, etc.):
# plt.title('Time Series of '+header['OBJECT'])
plt.grid(True)
plt.xlabel("Wavelength [Angstrom]")
# # Legend can be customized, e.g. font size (fontsize)
# plt.legend(
#     bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
#     loc=3,
#     ncol=4,
#     mode="expand",
#     borderaxespad=0.0,
#     fontsize=6,
# )

# plt.pause(20)  # uring this time, the graphic can be actively edited.


frage = input(
    'Would you like to save the graphic as png or pdf? Then please enter "png" or "pdf": ')
if frage == 'pdf':
    fig.savefig(files+'_Zeitserie.pdf')
if frage == 'png':
    fig.savefig(files+'_Zeitserie.png')


# Generate the spectrum section, plot and save as fit
min = 0
max = 10000
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1
    if header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"] >= min:
        min = header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"]
    if header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"] + \
            header["CDELT1"] * header["NAXIS1"] <= max:
        max = header["CRVAL1"] + \
            (1 - header["CRPIX1"]) * header["CDELT1"] +\
            header["CDELT1"] * header["NAXIS1"]
print("\nCommon wavelength range: ",
      int(min), ' bis ', int(max),  ' Angstrom')

print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Begin: "))
b = float(input("End: "))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    if header["CRVAL1"] <= a and header["CRVAL1"] +\
            header["CDELT1"] * header["NAXIS1"] >= b:
        aindex = int((a - header["CRVAL1"]) / header["CDELT1"])
        bindex = int((b - header["CRVAL1"]) / header["CDELT1"])
        filename = (
            filelist[i].rsplit(".")[0] + "_" + str(int(a)) +
            "_" + str(int(b)) + ".fits"
        )
        newflux = flux[aindex:bindex]
        header["CRVAL1"] = header["CRVAL1"] + aindex * header["CDELT1"]
        header["CRPIX1"] = 1
        header["NAXIS1"] = bindex - aindex
        print('Verwendet: ', i, header["CRVAL1"], header["NAXIS1"])
        fits.writeto(filename, newflux, header,
                     overwrite=True, output_verify="silentfix")

print('End of the program')
plt.close('all')
