#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads spectra catalog (time series of normalized 1d spectra in
tab format). Fits the selected line using nth degree regression
and determines the radial velocity RV from the minimum.
Plots all fittings and returns determined data as ascii files
(comma-separated, as .csv).

Stand 20221105

@author: lothar schanne
"""

import numpy as np
from astropy.io import ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pyplot as plt
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# local module, to be found in the "linelist" folder,
# must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print("Number of spectra: ", len(filelist), "\n")

# Selection of line or wavelength
linie, wellenlaenge = Linienlisten.linienauswahl()

frage_grafikenspeichern = input(
    'If you want to save the graphics, enter y: ')


# # Processing the filelist, reading flux and header, selecting the flux around the line:

# Definition of variables
RV = np.zeros(len(filelist))
apex = np.zeros(len(filelist))
linien_minwave = np.zeros(len(filelist))
miniwave = np.zeros(len(filelist))
miniflux = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="tab")

    # **********************************************************
    # Adjust the width of the search interval, width in angstroms.
    suchintervall = 3
    # ************************************************************

    intervall = np.extract(
        abs(ts.columns[0] - wellenlaenge) < suchintervall, ts)
    a = np.zeros(len(intervall))
    b = np.zeros(len(intervall))

    for j in range(len(intervall)):
        a[j], b[j] = intervall[j]
    # apex[i] = np.min(b)
    linien_argmin = np.argmin(b)
    # absolute line minimum (pixel by pixel)
    linien_minwave[i] = a[linien_argmin]

    fig = plt.figure(i)
    plt.plot(a, b, "-", linewidth=0.5)
    plt.title(filelist[i].rstrip(".dat") + "_" + linie)

    # Adjust the following interval [Angstroms], wide for noisy lines,
    # small for non-noisy, narrow lines.
    linienminimum_wave = a[abs(a - linien_minwave[i]) <= 0.5]
    linienminimum_flux = b[abs(a - linien_minwave[i]) <= 0.5]

    # Regression, adjust grade (2 or 4 or 6)
    grad = 4
    model = np.poly1d(np.polyfit(linienminimum_wave, linienminimum_flux, grad))

    polyline = np.linspace(linienminimum_wave[0], linienminimum_wave[-1], 100)
    modflux = model(polyline)

    # Calculate line minimum
    miniflux[i] = modflux.min()
    miniwave[i] = polyline[modflux.argmin()]

    # Plotting
    plt.plot(polyline, modflux, "--")
    plt.plot(miniwave[i], miniflux[i], "o", color="black")
    plt.pause(.1)
    if frage_grafikenspeichern == 'y':
        plt.savefig(filelist[i].rstrip(".dat") + "_" + linie + ".png")

    # RV without barycentric correction
    RV[i] = (miniwave[i] - wellenlaenge) / wellenlaenge * 299792

# # # Plot of RV's
# # fig = plt.figure()
# # plt.plot(obs_time, RV, 'bo', markersize=1)

# Saving as an ascii file
ascii.write(
    [filelist, RV, miniflux],
    linie + "_RV_Regression_grad" + str(grad) + ".csv",
    overwrite=True,
    names=["Spectrum               ",
           "        RV          ", "         miniFlux", ],
    format="csv",
)

fr = input('When you have finished viewing the graphics and want to end the\
 program, press the Enter key.')
plt.close('all')
print('\nProgram is completed')
