#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads a time series of 1d spectra in fits format.
A range of the Julian date can be selected, the spectra of which are plotted
plotted and a wavelength range.
Plot the spectra sections with an offset on top of each other and save the plot
as .png and .pdf.

Stand 20221105
@author: Lothar Schanne
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob

plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

# Create file list. All spectra in one (sub)folder.
files = input("Path and name of the files (use wildcards) : ")
filelist = glob.glob(files)

# Aalphabetical sorting. If the name is correct, this results in a
# temporal order
filelist.sort()

obs_time = np.zeros(len(filelist))
# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")
print('List of observation times of the spectrum series in JD: ')
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "JD-OBS" in header:
        obs_time[i] = float(header["JD-OBS"])
    elif 'JD' in header:
        obs_time[i] = float(header["JD"])
    if "JD-MID" in header:
        obs_time[i] = float(header["JD-MID"])
    elif "JD_OBS" in header:
        obs_time[i] = float(header["JD_OBS"])
    elif "MJD-OBS" in header:
        mjd = header["MJD-OBS"]
        obs_time[i] = mjd + 2400000.5
    elif "BAS_MJD" in header:
        obs_time[i] = float(header["BAS_MJD"]) + 2400000.5
    else:
        print("There is no observation time in the header of ",
              filelist[i])
    print(filelist[i], ": ", obs_time[i])

# Enter the range of observation times as JD
JD_anfang = float(input('Enter the start of the range of \
observation times (JD): '))
JD_ende = float(input('Enter the end of the range of \
observation times (JD): '))

# Input of the imaged wavelength range
lambda_anfang = float(input('Specify the beginning of the \
wavelength range to be imaged:'))
lambda_ende = float(input('Specify the end of the \
wavelength range to be imaged: '))

# Parameters for the distance between the spectra in the second plot
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")

fig = plt.figure(1)
zaehler = 0

# Read header and flux
for i in range(len(filelist)):
    print(filelist[i], ":")
    flux, header = fits.getdata(filelist[i], header=True)
    # Generate the spectrum sections:
    step = header["CDELT1"]
    refpix = header["CRPIX1"]
    # JD der Beobachtung:
    JD = obs_time[i]

    if JD >= JD_anfang and JD <= JD_ende:
        wave_erstesPix = header["CRVAL1"] - step * (refpix - 1)

        wave = np.zeros(header["NAXIS1"], dtype=float)
        wave_bereich = np.array([])
        flux_bereich = np.array([])
        for k in range(header["NAXIS1"]):
            wave[k] = wave_erstesPix + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])
        plt.plot(wave_bereich, flux_bereich + zaehler *
                 offset, "k-", label=str(JD), linewidth=1.)
        plt.text(wave_bereich[1], 1.0 + zaehler * offset, str(JD), ha="left",
                 size=7)
        zaehler = zaehler + 1

    else:
        pass


# Customize plot properties (grid, labels, etc.):
# plt.ylim(0.2, 1.0 + zaehler * offset + zusatz)
plt.title("Time Series of " + obj + ', JD = '+str(JD_anfang)+' bis '+str(JD_ende)
          )
plt.grid(True)
plt.xlabel("Wavelength in Angström")
fig.savefig("timeseries.pdf"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.pdf', format='pdf')

plt.pause(.1)

print('ZClick on the last opened diagram to exit the program.')
plt.waitforbuttonpress(-1)
plt.close('all')
