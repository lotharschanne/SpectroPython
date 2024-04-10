#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads spectra catalog (time series of normalized 1d spectra in
tab format). Determines the approximate line minimum by interaction and
determines the radial velocity RV from the minimum determined by regression.
If desired, a second line can be treated.
Plots all fittings and outputs determined data as ascii files (comma-separated,
as .csv).
The RV's are not barycentrically corrected.


Stand 20220408

@author: lothar schanne
"""

import numpy as np
from astropy.io import ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pyplot as plt

# local module, must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Print out list
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")

# Selection of the line or wavelength
linie, wellenlaenge = Linienlisten.linienauswahl()

frage_grafikenspeichern = input(
    'If you want to save the graphics, enter y: ')


# Processing the filelist, reading in flux and wavelengths,
# Selection of the flux around the line:

# Definition of variables
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
miniwave1 = np.zeros(len(filelist))
miniwave2 = np.zeros(len(filelist))
miniflux1 = np.zeros(len(filelist))
miniflux2 = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="tab")

    # **********************************************************
    # Adjust the width of the search interval, width in angstroms.
    suchintervall = 4
    # ************************************************************

    intervall = np.extract(
        abs(ts.columns[0] - wellenlaenge) < suchintervall, ts)
    a = np.zeros(len(intervall))
    b = np.zeros(len(intervall))
    for j in range(len(intervall)):
        a[j], b[j] = intervall[j]

    fig = plt.figure(i)
    plt.plot(a, b, "-", linewidth=0.5)
    plt.title(filelist[i].rstrip(".dat") + "_" + linie)

    print("\nSpectrum ", filelist[i])
    print("Click on the minimum of the first line: ")
    pts = np.asarray(plt.ginput(n=1, timeout=-1))

    # Adjust the following interval [Angstroms], wide for noisy lines,
    # small for non-noisy, narrow lines.
    linienminimum_wave1 = a[abs(a - pts[0, 0]) <= 0.5]
    linienminimum_flux1 = b[abs(a - pts[0, 0]) <= 0.5]

    ################ Regression, Adjust degree (2 oder 4 oder 6) ###########
    grad = 4
    ########################################################################
    model = np.poly1d(np.polyfit(
        linienminimum_wave1, linienminimum_flux1, grad))

    polyline = np.linspace(
        linienminimum_wave1[0], linienminimum_wave1[-1], 100)
    modflux = model(polyline)

    # Calculate line minimum
    miniflux1[i] = modflux.min()
    miniwave1[i] = polyline[modflux.argmin()]

    # Plotting
    plt.plot(polyline, modflux, "--")
    plt.plot(miniwave1[i], miniflux1[i], "o", color="black")
    # plt.show()
    # plt.savefig(filelist[i].rstrip(".dat") + "_" + linie + ".png")

    RV1[i] = (miniwave1[i] - wellenlaenge) / wellenlaenge * 299792

    frage = input(
        "Would you like to click on the minimum of a second line? Then enter y: "
    )
    if frage == "y":
        print("Click on the second line")
        pts = np.asarray(plt.ginput(n=1, timeout=-1))

        # Adjust the following interval [Angstroms], wide for noisy lines,
        # small for non-noisy, narrow lines.
        linienminimum_wave2 = a[abs(a - pts[0, 0]) <= 0.3]
        linienminimum_flux2 = b[abs(a - pts[0, 0]) <= 0.3]

        # Regression
        model = np.poly1d(np.polyfit(
            linienminimum_wave2, linienminimum_flux2, grad))

        polyline = np.linspace(
            linienminimum_wave2[0], linienminimum_wave2[-1], 100)
        modflux = model(polyline)

        # Calculate line minimum
        miniflux2[i] = modflux.min()
        miniwave2[i] = polyline[modflux.argmin()]

        # Plotting
        plt.plot(polyline, modflux, "--")
        plt.plot(miniwave2[i], miniflux2[i], "o", color="black")

        RV2[i] = (miniwave2[i] - wellenlaenge) / wellenlaenge * 299792
    else:
        RV2[i] = np.NaN

    plt.pause(.1)
    if frage_grafikenspeichern == 'y':
        plt.savefig(filelist[i].rstrip(".dat") + "_" + linie + ".png")


# Saving as an ascii file
ascii.write(
    [filelist, RV1, RV2],
    linie + "_RV_interactiv2lines" + ".csv",
    overwrite=True,
    names=["Spectrum", "RV1", "RV2"],
    format="csv",
)

fr = input('When you have finished viewing the graphics and want to end the\
 program, press the Enter key.')
plt.close('all')
print('\nProgram is finished')
