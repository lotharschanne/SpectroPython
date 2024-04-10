#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a spectra catalog (time series of normalized 1d spectra in tab format,
2 columns WAVE and FLUX). Determination of the line minimum via interaction and
determination of the radial velocity from the minimum.
Plots all spectra to mark up to 2 line minima (in the case of a double star
double star SB2) and outputs the determined data (RV and apex) as ascii files
(comma-separated, as .csv). The RV's are not barycentrically corrected.


Stand 20230610

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

# Printout of the list
print("\nSpectrum list: \n")
print("Number of spectra: ", len(filelist), "\n")

# Selection of the line or wavelength
linie, wellenlaenge = Linienlisten.linienauswahl()


# Processing the filelist, reading in flux and wavelengths, selecting the flux
# around the line:

# Definition of variables
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
apex1 = np.zeros(len(filelist))
apex2 = np.zeros(len(filelist))
linien_minwave1 = np.zeros(len(filelist))
miniwave1 = np.zeros(len(filelist))
miniwave2 = np.zeros(len(filelist))
miniflux1 = np.zeros(len(filelist))
miniflux2 = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="csv")

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
    RV1[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299792
    apex1[i] = pts[0, 1]
    print("First line, RV =", RV1[i], "Apex =", apex1[i])

    frage = input(
        "Would you like to click on the minimum of a second line? Then enter y: "
    )
    if frage == "y":
        print("Click on the second line")
        pts = np.asarray(plt.ginput(n=1, timeout=-1))
        RV2[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299792
        apex2[i] = pts[0, 1]
        print("Zweite Linie, RV =", RV2[i], "Apex =", apex2[i])
    else:
        RV2[i] = np.NaN
        apex2[i] = np.NaN

    plt.savefig(filelist[i].rstrip(".dat") + "_" + linie + ".png")
    plt.close()


# Saving as an ascii file
ascii.write(
    [filelist, RV1, apex1, RV2, apex2],
    linie + "_RV_interactive" + ".csv",
    overwrite=True,
    names=["Spectrum", "RV1", "Apex1", "RV2", "Apex2"],
    format="csv",
)

fr = input('When you have finished viewing the graphics and want to end the\
 program, press the Enter key.')
plt.close('all')
print('\nProgram is finished')
