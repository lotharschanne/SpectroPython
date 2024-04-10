#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates for a time series of spectra normalized to the continuum in the
fits format the equivalent widths of several contiguous areas of a line.
The integration limits are entered interactively.

The first spectrum is plotted and displayed for 15 seconds. During this
period, the graphics window is switched to interactive mode so that the spectrum
and the integration wavelength range for the EW calculation
can be selected optically. This reaction time window of 15 seconds can be changed in
line 75.

After selecting the wavelength range to be displayed for all spectra in the series
you will be asked for the number of line regions for which the EW is to be determined.
The individual spectra of the series are then displayed graphically and
the limits of the spectra to be evaluated separately are requested by clicking on the limits.
This allows the wavelengths of the integration limits are determined individually
for each spectrum in the series.
At the end, the determined EWs are written to an ascii file.

Version 20230524

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards): ")
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
print("\Spectrum list: \n")
print(filelist)
print("Number of spectra: ", len(filelist), "\n")

# Graphic first spectrum
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

print('Please edit the graphics window (zoom out the line \
enlarge) and then wait')
plt.pause(10)
print('Now select the wavelength range of the line with 2 mouse clicks for \
all spectra together')

pts = []
pts = np.asarray(plt.ginput(n=2, timeout=-1))
plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
bereichsanfang = pts[0, 0]
bereichsende= pts[-1, 0]

anzahl = 1 + int(input('Enter the desired number of line areas to be integrated separately: '))



EW = np.zeros([len(filelist), anzahl-1])
JD = np.zeros(len(filelist))

# Processing the filelist
for k in np.arange(len(filelist)):
    sp = fits.open(filelist[k], ignore_missing_end=True)

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

    WAVE = []
    FLUX = []

    for n in np.arange(len(wave)):
        if wave[n] >= bereichsanfang and wave[n] <= bereichsende:
            WAVE.append(wave[n])
            FLUX.append(flux[n])

    # Plot the spectrum
    fig = plt.figure()
    plt.plot(WAVE, FLUX)
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("ADU")
    plt.title("Spectrum " + filelist[k])
    plt.grid(True)


    # Interactive setting of the integration limits
    # by clicking on the spectrum from left to right
    pts = []
    pts = np.asarray(plt.ginput(n=anzahl, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[-1, 0]
    print("Selected integration area:", begin, " bis", end)

    for m in range(len(pts)-1):
        for n in np.arange(len(WAVE)):
            if WAVE[n] <= pts[m,0] and WAVE[n + 1] > pts[m,0]:
                begin_n = n
                break
        for n in np.arange(len(WAVE)):
            if WAVE[n] <= pts[m+1,0] and WAVE[n + 1] > pts[m+1,0]:
                end_n = n
                break
        ew = 0
        for o in np.arange(end_n - begin_n):
            dif = 1 - FLUX[begin_n + o]
            ew += dif
        EW[k, m] = ew * sp[0].header["CDELT1"]

    print(filelist[k], "   EW =", EW[k])
    print()


ascii.write(
    [filelist, JD, EW],
    'EW_' + str(int(begin)) + "_" + str(int(end)) + '.ecsv',
    names=["Spectrum", 'JD', "EW"],
    overwrite=True,
    format="ecsv",
)



print('To exit the program, click on the last opened diagram')
plt.waitforbuttonpress(-1)
plt.close('all')
