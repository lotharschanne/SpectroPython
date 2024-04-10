#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a spectrum catalog (time series of normalized 1d spectra in the
fits format). Calculates the heliocentric correction from the observation time
and the (to be adjusted) coordinates of the observer and object,
fits the selected line in the minimum range using regression and uses the
heliocentrically corrected minimum to determine the heliocentrically corrected
radial velocity RV.
Plots and saves all fittings and outputs determined data as ascii files
(comma-separated, as .csv).

Stand 20221105

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt
from PyAstronomy import pyasl
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# local module, to be found in the "Blocks" folder,
# must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


# ************* Coordinates of the observatory in the following form: *************
# longitude = 289.5967661 in degrees oder [grad, min, sec]
# latitude = -24.62586583 in degrees oder [grad, min, sec]
# altitude = 2635.43 in Meter

# ***************  PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# If this is omitted, the barycentric corrections are wrong

# Coordinates of Berthold
# longitude = +7.4775
# latitude = 49.47527777778
# altitude = 200
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]


# Coordinates of Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Coordinates of Siegfried Hold
longitude = +15.68461111
latitude = 47.00161111
altitude = 380


# Important !!!!!!!!!!!!!!!! :
# esign : int, optional, {-1,0,1}
# Explicit sign with -1 representing negative sign, +1 representing positive
# sign, and 0 indicating no explicit sign specification. The explicit sign is
# necessary if negative southern coordinates are specified but d is 0 and,
# thus, cannot carry the sign.

if type(longitude) == list:
    longitude = pyasl.dmsToDeg(longitude[0], longitude[1], longitude[2])

if type(latitude) == list:
    latitude = pyasl.dmsToDeg(latitude[0], latitude[1], latitude[2], esign=+1)

# ********************************************************************

# ***************  PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# ************ Enter the coordinates of the star in the following form: *******
# If this is omitted, the barycentric corrections are incorrect
# ra2000 = 030.20313477 in degrees
# dec2000 = -12.87498346 in degrees
# or as coordinate string RA DEC in the form "hh mmm ss +dd mm ss"

# # Coordinates of del Cep
# ra2000 = 337.29277083333335
# dec2000 = +58.415198
# coord = "22 29 10.265 +58 24 54.714"

# Coordinates of gam Cyg
# ra2000 = 305.55708
# dec2000 = +40.2566

# Coordinates of Betageuze
# ra2000 = 88.7929583
# dec2000 = 7.40705555

# Coordinates of theta1 Ori C
# ra2000 = 83.81858333
# dec2000 = -5.389694444

# Coordinates of 7 And
# ra2000 = 348.1375
# dec2000 = 49.40620275
# coord = "23 12 33 +49 24 22.3299"

# Coordinates of gam Cyg
# ra2000 = 305.55708
# dec2000 = 40.25666666

# Coordinates of OX Aurigae
# ra2000 = 103.25
# dec2000 = 38.86916
# coord = "06 53 01.41099 +38 52 08.9353"

# Coordinates of Polaris (alp UMi)
# coord = "02 31 49 +89 15 51"

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Reading in the star coordinates via the Internet
Frage = input(
    'Would you like to search for the star coordinates on the Internet? Then enter "y": ')
if Frage == 'y':
    star = input('Enter the name of the object star: ')
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************

# Create file list. Spectra in a (sub)folder.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print("Number of spectra: ", len(filelist), "\n")

# Selection of the line or wavelength
linie, wellenlaenge = Linienlisten.linienauswahl()

# Entering the system speed
systemgeschwindigkeit = float(
    input("Enter a system velocity in km/s: ")
)
systemgeschwindigkeit = systemgeschwindigkeit / 299772 * wellenlaenge

frage_bary = input(
    'Would you like to correct the RV barycentrically? Then enter "y"')

frage_grafikenspeichern = input(
    'If you want to save the graphics, enter y: ')


# Processing the filelist, reading flux and header, selecting the flux around the line:

# Definition of variables
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
miniflux = np.zeros(len(filelist))

#################################################################
# Regression, adjust level (2 or 4 or 6)
grad = 2
#################################################################

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    if "JD-OBS" in header:
        obs_time[i] = float(header["JD-OBS"])
    elif 'JD' in header:
        obs_time[i] = float(header["JD"])
    elif "JD_OBS" in header:
        obs_time[i] = float(header["JD_OBS"])
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = float(header["BAS_MJD"])
    else:
        print("There is no observation time in the header of ",
              filelist[i])
        break
    # Calculation of the heliocentric correction and time:
    if frage_bary == 'y':
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i], debug=False
        )
    else:
        corr[i] = 0

    print("\n" + filelist[i] + ":")
    print("Date of observation: ", obs_time[i])
    print("Barycentric correction [km/s]: ", corr[i].round(2))

    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1

    step = float(header["CDELT1"])

    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step

    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

    # **********************************************************
    # Adjust the width of the search interval, width in angstroms.
    suchintervall = int(6 / step)
    # ************************************************************
    intervallflux = np.zeros(suchintervall)
    intervallwave = np.zeros(suchintervall)
    for j in range(suchintervall):
        intervallflux[j] = flux[j +
                                index_wellenlaenge - int(suchintervall / 2)]
        intervallwave[j] = (
            lambda0 + (j + index_wellenlaenge - int(suchintervall / 2)) * step
        )

    linienminimum_wave_index = intervallflux.argmin()

    ##########################################################################
    # Adjust the following interval, wide for noisy lines (-20,+21)
    # small for non-noisy, narrow lines (-2,+3)
    ##########################################################################
    linienminimum_intervall = np.arange(
        linienminimum_wave_index - 5, linienminimum_wave_index + 6
    )

    # Wavelengths and flux directly around the line minimum
    wl = intervallwave[linienminimum_intervall]
    fl = intervallflux[linienminimum_intervall]

    model = np.poly1d(np.polyfit(wl, fl, grad))

    polyline = np.linspace(wl[0], wl[-1], 100)
    modflux = model(polyline)

    # Calculate line minimum
    miniflux[i] = modflux.min()
    miniwave = polyline[modflux.argmin()]
    miniwave_bc = miniwave * (1 + corr[i] / 299792)

    # Plotting
    fig = plt.figure()
    plt.plot(intervallwave, intervallflux, 'o')
    plt.plot(polyline, modflux)
    plt.plot(miniwave, miniflux[i], "o", color="black")
    plt.title(filelist[i])
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("relative Intensity")
    plt.pause(2)
    if frage_grafikenspeichern == 'y':
        plt.savefig(filelist[i].rstrip(".fit") +
                    "_" + linie + "_Regression.png")
    plt.close()

    # RV without barycentric correction
    RV[i] = (miniwave - wellenlaenge) / wellenlaenge * 299792
    print()
    print(filelist[i]+" nicht korrigierte RV: ", RV[i].round(2))
    # RV only bc-corrected, sysv not taken into account:
    if frage_bary == 'j':
        RV_bc[i] = (miniwave_bc - wellenlaenge) / wellenlaenge * 299792
        print(filelist[i]+" barycentrically corrected RV: ", RV_bc[i].round(2))
    else:
        RV_bc[i] = None


# # Plot of RV's
# fig = plt.figure()
# plt.plot(obs_time, RV, 'bo', markersize=1)

# Saving as an ascii file
ascii.write(
    [filelist, obs_time, corr, RV, RV_bc, miniflux],
    linie + "_RV_Regression_grad" + str(grad) + ".csv",
    overwrite=True,
    names=[
        "Spectrum",
        "JD",
        "BC",
        "RV",
        "RV_bc",
        "minimumFlux",
    ],
    format="csv",
)

plt.close('all')
print('\nThe program is finished')
