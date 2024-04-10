#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads spectra catalog (time series of normalized 1d spectra in
fits format). Calculates the heliocentric correction from the time of
observation and the coordinates of the observer and object, fits the selected line per spline and uses the heliocentrically corrected minimum to determine the
heliocentrically corrected radial velocity RV from the heliocentrically
corrected minimum.
Plots and saves all fittings and outputs determined data as ascii files
(comma-separated, as .csv).

Stand 20220408

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from PyAstronomy import pyasl
from scipy import interpolate
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# local module, to be found in the folder "lienelist",
# must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# ************* Coordinates of the observatory in the following form: *************
# longitude = 289.5967661 in Grad or [grad, min, sec]
# latitude = -24.62586583 in Grad or [grad, min, sec]
# altitude = 2635.43 in Meter

# ***************  PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# If this is omitted, the barycentric corrections are wrong

# Coordinates of Berthold
# longitude = +7.4775
# latitude = 49.47527777778
# altitude = 200
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]


# Coordinates of  Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Coordinates of  Siegfried Hold
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

# *************** PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# ************ Enter the coordinates of the star in the following form: *******
# If this is omitted, the barycentric corrections are incorrec
# ra2000 = 030.20313477 in Grad
# dec2000 = -12.87498346 in Grad
# or as a coordinate string RA DEC in the form "hh mmm ss +dd mm ss"

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
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

#Print out the list
print("\nSpectra list: \n")
print("Number of spectra: ", len(filelist), "\n")

# Auswahl der Linie oder Wellenlänge
linie, wellenlaenge = Linienlisten.linienauswahl()

frage_bary = input(
    'Would you like to correct the RV barycentrically? Then enter "y"')

frage_grafikenspeichern = input(
    'If you want to save the graphics, enter y: ')

# Processing the filelist, reading flux and header, selecting the flux around the line:

# Definition of variables
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
linienminimum_flux = np.zeros(len(filelist))
linienminimum_wave = np.zeros(len(filelist))
linienminimum_bc = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))


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
        hjd[i] = 0

    print("\n" + filelist[i] + ":")
    print("Date of observation: ", obs_time[i])
    print("Barycentric correction [km/s]: ", corr[i])

    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1

    step = float(header["CDELT1"])

    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step

    index_wellenlaenge = int((wellenlaenge - lambda0) / step)

    wave = np.zeros(header['NAXIS1'], dtype=float)
    for k in range(header['NAXIS1']):
        wave[k] = lambda0 + k * step

    # **********************************************************
    # Adjust the width of the search interval, width in angstroms.
    suchintervall = int(8 / step)
    # ************************************************************
    intervallflux = np.zeros(suchintervall)
    intervallwave = np.zeros(suchintervall)
    for j in range(suchintervall):
        intervallflux[j] = flux[j +
                                index_wellenlaenge - int(suchintervall / 2)]
        # intervallwave[j] = (
        #     lambda0 + (j + index_wellenlaenge - int(suchintervall / 2)) * step
        intervallwave[j] = wave[j +
                                index_wellenlaenge - int(suchintervall / 2)]

    # Spline over the spectrum in the interval
    tck = interpolate.splrep(intervallwave, intervallflux, s=0.0001)
    newwaveintervall = np.arange(
        intervallwave[0], intervallwave[-1], step / 50)
    newfluxinterpolated = interpolate.splev(newwaveintervall, tck, der=0)
    linienminimum_flux[i] = newfluxinterpolated.min()
    # linienminimum_wave[i] = intervallwave[0] + \
    #     newfluxinterpolated.argmin()*step/10
    linienminimum_wave[i] = newwaveintervall[newfluxinterpolated.argmin()]
    linienminimum_bc[i] = linienminimum_wave[i] * (1 + corr[i] / 299792)

    RV[i] = (linienminimum_wave[i] - wellenlaenge) / wellenlaenge * 299792
    print('\n' + filelist[i])
    print('uncorrected RV: ', RV[i])
    # RV only bc-corrected, sysv not taken into account:
    RV_bc[i] = (linienminimum_bc[i] - wellenlaenge) / wellenlaenge * 299792
    print("barycentrically corrected RV: ", RV_bc[i])

    # Plotten
    fig = plt.figure()
    plt.title(filelist[i] + ', JD: ' + str(obs_time[i]), fontsize=7)
    plt.plot(intervallwave, intervallflux, "o")
    plt.plot(newwaveintervall, newfluxinterpolated, "-")
    plt.plot(linienminimum_wave[i], linienminimum_flux[i], "o", color="red")
    plt.pause(1)
    if frage_grafikenspeichern == 'y':
        plt.savefig(filelist[i].rstrip(".fit") + "_" + linie + "_Spline.png")
    plt.close()

# # Plot of RV's
# fig = plt.figure()
# plt.plot(obs_time, RV, 'bo', markersize=1)

# Saving as an ascii file
ascii.write(
    [filelist, obs_time, corr, RV, RV_bc, linienminimum_flux],
    linie + "_RV_perSpline" + ".csv",
    overwrite=True,
    names=[
        "      Spectrum     ",
        "      JD      ",
        "     BC     ",
        "     RV     ",
        "      RV_bc     ",
        "      miniFlux",
    ],
    format="csv",
)


# ascii.write([obs_time, apex3], linie+'_apex'+'.dat', overwrite=True,
#             names=['JD', 'apex'], format='tab')


# # Saving obs_time and corr in ascii-file
# ascii.write([obs_time, corr], linie+'_table_obs_time_bc.dat', overwrite=True,
#             names=['JD', 'BARYCORR'], format='tab')

# # Saving obs_time and jhd in ascii-file
# ascii.write([obs_time, hjd], linie+'_table_obs_time_hjd.dat', overwrite=True,
#             names=['JD', 'HJD'], format='tab')


plt.close('all')
print('\nThe program is finished')
