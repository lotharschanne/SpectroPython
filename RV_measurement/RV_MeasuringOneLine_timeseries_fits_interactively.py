#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads spectra catalog (time series in fits format).
Calculates from the observation time and the (to be adjusted) coordinates of
the observer and object, calculates the heliocentric correction, shows the
selected line in a plot, the minimum of the selected line is to be clicked and
determines the heliocentrically corrected radial velocity RV_bc.
Outputs determined data as csv-ascii files (comma-separated).

Stand 20221105
@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt
from PyAstronomy import funcFit as fuf
from PyAstronomy import pyTiming as pyt
from PyAstronomy import pyasl
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from time import sleep

# local module, must be in the same directory as the script (or the
# directory of the module in the Python path so that Python can find it)
import Linienlisten


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend


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


# Coordinates of  Wise Observatory in Israel
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
# ra2000 = 030.20313477 in Grad
# dec2000 = -12.87498346 in Grad
# or as a coordinate string RA DEC in the form "hh mmm ss +dd mm ss"

# # Coordinates ofdel Cep
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

# Ausdruck der Liste
print("\nSpectrum list: \n")
print("Number of spectra: ", len(filelist), "\n")


# Selection of the line or wavelength
linie, laborwellenlaenge = Linienlisten.linienauswahl()

# Enter the system velocity
sys = input(
    "Would you like to enter a system speed? Enter the number in km/s or n: "
)
if sys == "n":
    sysRVlambda = 0.0
else:
    sysRVlambda = float(sys) * laborwellenlaenge / 299792

# Processing the filelist, reading flux and header, selecting the flux around the line:
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
apex = np.zeros(len(filelist))

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "JD" in header:
        obs_time[i] = header["JD"]
    elif "JD-OBS" in header:
        obs_time[i] = header["JD-OBS"]
    elif "JD_OBS" in header:
        obs_time[i] = header["JD_OBS"]
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = header["BAS_MJD"]
    else:
        print("!!!!! No observation time in the header  !!!!")

    # Calculation of the heliocentric correction:
    corr[i], hjd[i] = pyasl.helcorr(
        longitude, latitude, altitude, ra2000, dec2000, obs_time[i], debug=False
    )

    print("\n" + filelist[i] + ":")
    print("Date of observation: ", obs_time[i])
    print("Barycentric correction [km/s]: ", corr[i])
    # print("Heliocentric Julian day: ", hjd[i])
    lambda0 = header["CRVAL1"]
    if "CRPIX1" in header:
        refpix = header["CRPIX1"]
    else:
        refpix = 1
    step = header["CDELT1"]

    wave = np.ones(header["NAXIS1"], dtype=float)
    header["CRVAL1"] = header["CRVAL1"] + (1 - refpix) * step
    for j in range(header["NAXIS1"]):
        wave[j] = header["CRVAL1"] + j * step + sysRVlambda

    # Plot the spectrum
    fig = plt.figure()
    plt.plot(wave, flux)
    plt.xlabel("Wavelength [Angström]", fontsize="medium")
    plt.ylabel("Flux")
    # plt.ylim(0.7, 1.03)
    plt.title("Spectrum " + filelist[i])
    plt.grid(True)
    plt.xticks(fontsize="small")
    plt.yticks(fontsize="small")
    plt.xlim(laborwellenlaenge + sysRVlambda - 5,
             laborwellenlaenge + sysRVlambda + 5)

    print("Click on the minimum of the line and then on the base (the continuum): ")
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))

    RV[i] = (pts[0, 0] - laborwellenlaenge) / laborwellenlaenge * 299792
    RV_bc[i] = RV[i] + corr[i]
    apex[i] = pts[0, 1] / pts[1, 1]
    print("\nBarycentric (and systemic) corrected RV:",
          RV_bc[i], "Apex:", apex[i])
    plt.pause(.1)


# # Plot of RV's
# fig=plt.figure()
# plt.plot(obs_time, RV,'bo', markersize=1)
# plt.plot(obs_time, RV_bc,'r+', markersize=1)

# Saving as an ascii file
ascii.write(
    [filelist, obs_time, RV_bc],
    linie + "_RV_bc" + ".csv",
    overwrite=True,
    names=['Spectrum',  "JD",  "RV_bc"],
    format="csv",
)
ascii.write(
    [filelist, obs_time, apex],
    linie + "_apex" + ".csv",
    overwrite=True,
    names=['Spectrum', "JD", "apex"],
    format="csv",
)
# Speichern von obs_time und corr in ascii-file
ascii.write(
    [filelist, obs_time, corr],
    linie + "_Tabelle_obs_time_bc.csv",
    overwrite=True,
    names=['Spectrum', "JD", "BARYCORR"],
    format="csv",
)

fr = input('When you have finished viewing the graphics and want to end the\
 program, press the Enter key.')
plt.close('all')
print('\nThe program is finished')
