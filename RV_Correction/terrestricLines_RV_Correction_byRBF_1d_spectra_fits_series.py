#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Input: Read in a spectrum series in fits format. These must not
be heliocentrically corrected.

The script corrects the calibration by determining the line minimum of known
terrestrial lines and converts the spectra with the determined mean RV if it
exceeds 3 km/s in magnitude.
Save the corrected spectra as fits. Showing the modeled terrestrial lines as a
plot. Print out the determined RVs. Writing of the RVs, the mean
value and the standard deviation in an ascii file.

@author: Lothar Schanne
Stand 20221116
"""

import numpy as np
from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pylab as plt
from scipy.interpolate import Rbf


# List of terrestrial lines used:
Linien = [
    6543.912,
    6552.646,
    6574.88]
print("List of terrestrial lines used:: ",
      Linien)

object = input('Enter the name of the object (without blanks): ')

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
miniwave = np.zeros([len(filelist), len(Linien)])
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
        print("number of values (abscissa), NAXIS1:     ", nax)
    else:
        print("NAXIS1 missing in header !")
    if "CRVAL1" in header:
        crval = header["CRVAL1"]
        print("Starting wavelength, CRVAL1:             ", crval)
    else:
        print("CRVAL1 missing in header !")
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
        suchintervall = int(1 / cdel)
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

        # Radial basis function (RBF) over the spectrum in the interval
        # adjust smooth, 0. = function goes through all points, >0. = equalize
        rbf = Rbf(intervallwave, intervallflux, smooth=0.1)

        newwaveintervall = np.arange(
            intervallwave[0], intervallwave[-1], cdel / 10)
        newfluxinterpolated = rbf(newwaveintervall)
        miniflux[i, m] = newfluxinterpolated.min()
        # linienminimum_wave[i] = intervallwave[0] + \
        #     newfluxinterpolated.argmin()*step/10
        miniwave[i, m] = newwaveintervall[newfluxinterpolated.argmin()]

        RV[i, m] = (miniwave[i, m] - wellenlaenge) / wellenlaenge * 299792

        # Plotten
        fig = plt.figure()
        plt.plot(intervallwave, intervallflux, "o")
        plt.plot(newwaveintervall, newfluxinterpolated, "-")
        plt.plot(miniwave[i, m], miniflux[i, m], "or")
        # plt.savefig(filelist[i].rstrip(".fit") + "_" + Linien[m] + "_RBF.png")

        # RV
        RV[i, m] = (miniwave[i, m] - wellenlaenge) / wellenlaenge * 299792
        print("RV: ", RV[i, m].round(2))

    RV_mean[i] = RV[i].mean()
    RV_std[i] = RV[i].std()
    print('Mean RV: ', RV_mean[i], '\nStandard deviation RV: ', RV_std[i])

    # Shift the spectrum under the condition that |RV_mean| > 3 km/s
    if abs(RV_mean[i]) >= 3.:
        flux_rv, wave_rv = pyasl.dopplerShift(
            wave, flux, -RV_mean[i], edgeHandling="firstlast")

        # Save the RV-corrected spectrum in ascii-file
        # ascii.write(
        #     [wave, flux_rv],
        #     filelist[i].rstrip(".fits") + "_perRBF_RVcorrected.dat",
        #     overwrite=True,
        #     names=["WAVE", "FLUX"],
        #     format="tab",
        # )

        # Writing the RV-corrected spectrum to fits-file
        header["CRVAL1"] = wave[0]
        header["CRPIX1"] = 1
        header["NAXIS1"] = len(wave)
        newfile = filelist[i].rstrip(".fits") + "_perRBF_RVcorrected.fits"
        fits.writeto(newfile, flux_rv, header, overwrite=True,
                     output_verify="silentfix")

        print("New starting wavelength CRVAL1: ", header["CRVAL1"], "\n\n")

print('\nMean RVs of the spectra from the terrestrial lines: ', RV_mean, '\n')

# Save the RVs as an ascii file (csv):
ascii.write(
    [filelist, RV[:, 0], RV[:, 1], RV[:, 2], RV_mean, RV_std],
    object + "_terrLinesCorrigated" + ".csv",
    overwrite=True,
    names=[
        "Spectrum",
        "RV[0]", 'RV[1]', 'RV[2]',
        "RV_mean", 'RV_std'
    ],
    format="csv",
)
