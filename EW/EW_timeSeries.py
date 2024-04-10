#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates for a time series on the continuum of normalized spectra in the
fits format the equivalent width EW of a line.
Input of the integration limits graphically-interactively or manually.
The EW calculation assumes that all spectra in the series have the same
wavelength interval can be used for the calculation of the integral,
no significant RV changes take place, or if they do, that the line is isolated
(i.e. the flux = 1 in its environment). It is best to use barycentrically
corrected spectra.
The first spectrum is plotted and displayed for 15 seconds. During this
period, the graphics window is interactive, so that you can enlarge the spectrum
and the integration wavelength range for the EW calculation can be selected optically.
This reaction time of 20 seconds can be changed in line 66.

Version 20230524

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\List of spectra: \n")
print(filelist)
print("Number of spectra: ", len(filelist), "\n")

# Plot first spectrum
sp = fits.open(filelist[0], ignore_missing_end=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

for i in np.arange(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
    )
    # The list wave contains the wavelengths of the pixels.
# Close the fits-file:
sp.close()

# Plot the spectrum
fig = plt.figure()
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("ADU")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
print('Please wait or edit the graphics window (e.g. enlarge a line)')
plt.pause(15)
plt.show(block=False)


frage1 = input(
    "Would you like to enter the integration limits numerically (m) or \
by mouse click (graphically, g)? Enter m or g: "
)

if frage1 == "g":
    # Interactive definition of the integration limits
    # If the display of the spectrum is interactively enlarged, delete the input points
    # with the right mouse button until the left side of the line to be
    # measured is clicked.
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    print("Selected integration range:", begin, " to", end)
    print()

if frage1 == "m":
    # Entering the integration limits
    begin = float(
        input("Enter the short-wave wavelength integration limit: ")
    )
    end = float(
        input("Enter the long-wave wavelength integration limit: "))
    print()

EW = np.zeros(len(filelist))
JD = np.zeros(len(filelist))

# Processing the filelist
for k in np.arange(len(filelist)):
    sp = fits.open(filelist[k], ignore_missing_end=True)
    # Header verification
    if "JD-OBS" in sp[0].header:
        JD[k] = float(sp[0].header["JD-OBS"])
    elif 'JD' in sp[0].header:
        JD[k] = float(sp[0].header["JD"])
    elif "JD_OBS" in sp[0].header:
        JD[k] = float(sp[0].header["JD_OBS"])
    elif "BAS_MJD" in sp[0].header:
        JD[k] = float(sp[0].header["BAS_MJD"])
    else:
        print("No observation time in the header")
        JD[k] = 0

    try:
        sp[0].header["CRPIX1"]
    except:
        sp[0].header["CRPIX1"] = 1

    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )

    # Close the fits-file
    sp.close()

    for n in np.arange(len(wave)):
        if wave[n] <= begin and wave[n + 1] > begin:
            begin_n = n
            begin_wave = wave[n]
            begin_flux = np.median(flux[n - 2: n + 2])
            break
    for n in np.arange(len(wave)):
        if wave[n] <= end and wave[n + 1] > end:
            end_n = n
            end_wave = wave[n]
            end_flux = np.median(flux[n - 2: n + 2])
            break
    faktor = (end_flux - begin_flux) / (end_n - begin_n)
    ew = 0
    for m in np.arange((end_n - begin_n)):
        dif = begin_flux + m * faktor - flux[begin_n + m]
        ew += dif
    EW[k] = ew * sp[0].header["CDELT1"]

    print(filelist[k], "   EW =", EW[k])


ascii.write(
    [filelist, JD, EW],
    'EW_' + str(int(begin)) + "_" + str(int(end)) + '.dat',
    names=["Spectrum", 'JD', "EW"],
    overwrite=True,
    format="tab",
)

fig = plt.figure()
plt.stem(JD, -EW)
plt.xlabel("JD")
plt.ylabel("-EW in Angstroem")
plt.title('EW ' + str(begin) + ' to ' + str(end) + ' Angstroem')
plt.grid(True)
plt.savefig('EW_' + str(begin) + '_' + str(end) + '.png', format='png')
# plt.savefig('EW_' + str(begin) + '_' + str(end) + '.pdf', format='pdf')

print('To exit the program, click on the last opened diagram')
plt.waitforbuttonpress(-1)
plt.close('all')
