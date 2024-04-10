#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: Read in a spectrum series in fits format. These must not
be heliocentrically corrected.

The script corrects the calibration by determining the line minimum of known
terrestrial lines by regression and recalculates the spectra with the
determined mean RV if it exceeds 3 km/s in terms of magnitude.
Saving the corrected spectra as fits. Showing the modeled terrestrial
lines as a graphic. Print out the calculated RVs. Writing
the RV's, the mean value and the standard deviation into an ascii-file.

@author: Lothar Schanne
Stand 20221117
"""

import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pylab as plt

#################################################################
# Regression, adjust Grad (2 oder 4 oder 6)
grad = 2
#################################################################

# List of selectable terrestrial lines:
Linien = [
    6543.912,
    6552.646,
    6574.88]
print("List of selectable terrestrial lines:: ",
      Linien)

object = input('Enter the name of the object: ')

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")


miniflux = np.zeros([len(filelist), len(Linien)])
RV = np.zeros([len(filelist), len(Linien)])
RV_mean = np.zeros(len(filelist))
RV_std = np.zeros(len(filelist))


# Processing the filelist:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    print("\n", filelist[i])
    print("Original begin of wavelength scale, CRVAL1: ", header["CRVAL1"])

    if "NAXIS" in header:
        print("Dimension, NAXIS:                        ", header["NAXIS"])
    else:
        print("This is not a 1d spectrum !")
    if "NAXIS1" in header:
        nax = header["NAXIS1"]
        print("Number of values (abscissa), NAXIS1:     ", nax)
    else:
        print("NAXIS1 missing in header !")
    if "CRVAL1" in header:
        crval = header["CRVAL1"]
        print("Starting wavelength, CRVAL1:             ", crval)
    else:
        print("CRVAL1 missing in header!")
    if "CRPIX1" in header:
        crpix = header["CRPIX1"]
        print("Reference pixel, CRPIX1:             ", crpix)
    else:
        print("CRPIX1 missing in header !")
        crpix = 1
    if "CDELT1" in header:
        cdel = header["CDELT1"]
        print("Step size of the wavelength, CDELT1:    ", cdel)
    else:
        print("CDELT1 missing in header !")

    #   Creation of numpy arrays with the wavelengths and fluxes of the spectrum
    wave = np.ones(nax, dtype=float)
    for k in range(nax):
        wave[k] = crval + (k - crpix + 1) * cdel

    # Calculating the mean calibration error using the terrestrial lines

    for m in range(len(Linien)):
        wellenlaenge = Linien[m]
        index_wellenlaenge = int((wellenlaenge - wave[0]) / cdel)
        # **********************************************************
        # Adjust the width of the search interval, width in angstroms.
        suchintervall = int(3 / cdel)
        # ************************************************************
        intervallflux = np.zeros(suchintervall)
        intervallwave = np.zeros(suchintervall)
        for j in range(suchintervall):
            intervallflux[j] = flux[j +
                                    index_wellenlaenge - int(suchintervall / 2)]
            intervallwave[j] = (
                wave[0] + (j + index_wellenlaenge -
                           int(suchintervall / 2)) * cdel
            )

        linienminimum_wave_index = intervallflux.argmin()

        ##########################################################################
        # Adjust the following interval, wide for noisy lines (-20,+21)
        # small for non-noisy, narrow lines (-2,+3)
        ##########################################################################
        linienminimum_intervall = np.arange(
            linienminimum_wave_index - 1, linienminimum_wave_index + 2
        )

        # Wavelengths and flux directly around the line minimum
        wl = intervallwave[linienminimum_intervall]
        fl = intervallflux[linienminimum_intervall]

        model = np.poly1d(np.polyfit(wl, fl, grad))

        polyline = np.linspace(wl[0], wl[-1], 100)
        modflux = model(polyline)

        # Calculate line minimum
        miniflux[i, m] = modflux.min()
        miniwave = polyline[modflux.argmin()]
        print('')

        # Plotting
        fig = plt.figure()
        plt.plot(intervallwave, intervallflux, "o-")
        plt.plot(polyline, modflux)
        plt.plot(miniwave, miniflux[i, m], "o", color="black")
        plt.title(filelist[i] + 'Linie' + str(Linien[m]))
        plt.xlabel("Wavelength [Angström]")
        plt.ylabel("relative Intensity")
        # plt.savefig(filelist[i].rstrip(".fit") + "_" + Linie + "_Regression.png")

        # RV
        RV[i, m] = (miniwave - wellenlaenge) / wellenlaenge * 299792
        print("RV: ", RV[i, m].round(2))

    RV_mean[i] = RV[i].mean()
    RV_std[i] = RV[i].std()

    # # Shift that spectrum
    # flux_rv, wave_rv = pyasl.dopplerShift(
    #     wave, flux, -RV_mean, edgeHandling="firstlast")

    # ascii.write(
    #     [wave, flux_rv],
    #     filelist[i].rstrip(".fit") + "_RVcorrected.dat",
    #     overwrite=True,
    #     names=["WAVE", "FLUX"],
    #     format="tab",
    # )
    # # Writing the RV-corrected spectrum to fits-file
    # header["CRVAL1"] = wave[0]
    # header["CRPIX1"] = 1
    # header["NAXIS1"] = len(wave)
    # newfile = filelist[i].rstrip(".fit") + "_RVcorrected.fit"
    # fits.writeto(newfile, flux_rv, header, overwrite=True,
    #              output_verify="silentfix")

    # print("New starting wavelength CRVAL1: ", header["CRVAL1"], "\n\n")
print('\nMean RVs of the spectra from the terrestrial lines: \n', RV_mean, '\n')

# Save the RVs as an ascii file (csv):
ascii.write(
    [filelist, RV[:, 0], RV[:, 1], RV[:, 2], RV_mean, RV_std],
    object + "_terrLines_corrigated_byRegression" + ".csv",
    overwrite=True,
    names=[
        "Spectrum",
        "RV[0]", 'RV[1]', 'RV[2]',
        "RV_mean", 'RV_std'
    ],
    format="csv",
)

print('To exit the program, click on the last opened diagram.')
plt.waitforbuttonpress(-1)
plt.close('all')
