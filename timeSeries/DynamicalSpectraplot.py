#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generates a colored dynamic spectrum plot from a sequence of (barycentrically
corrected) fits-1d spectra, where the ordinates are the observation times
(a header entry named 'JD' must exist in the header of each spectrum!)
The wavelength range shown is selected.


20231212
@author: lothar
"""

import numpy as np
from astropy.io import fits
import glob
import matplotlib.pylab as plt


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the fits spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Ausdruck der Liste
print("\List of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

jd = np.zeros(len(filelist))
fluxmin = np.zeros(len(filelist))
fluxmax = np.zeros(len(filelist))
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    try:
        jd[i] = header['JD']
    except:
        print(filelist[i], ': No JD in the header')
        pass
print('JD range of the spectra extends from ', jd.min(), ' to ', jd.max())

# Enter the range of observation times as JD
JD_anfang = float(input('Enter the start of the range of \
observation points (JD) that you want to map: '))
JD_ende = float(input('Enter the end of the range of \
observation times (JD): '))

# Enter the wavelength range to be mapped
lambda_anfang = float(input('Specify the beginning of the \
wavelength range to be imaged: '))
lambda_ende = float(input('Specify the end of the \
wavelength range to be imaged: '))

obj = input("Please enter the object name: ")

plt.figure(1, figsize=(7, 10))
plt.title("Line development " + obj)
plt.grid(True)
plt.xlabel('Wavelength [Angstrom]')
plt.ylabel('JD')

# Determination of the minimum and maximum fluxes
flux_bereich_max = np.zeros(len(filelist))
flux_bereich_min = np.zeros(len(filelist))

for i in range(len(filelist)):
    print(i, filelist[i], ' is being examined')
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']
    JD_float = float(JD)
    wave_bereich = np.array([])
    flux_bereich = np.array([])

    if (JD_float >= JD_anfang) and (JD_float <= JD_ende):
        wave = np.zeros(header['NAXIS1'])
        for k in range(header['NAXIS1']):
            wave[k] = lambda0 + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])
        flux_bereich_max[i] = flux_bereich.max()
        flux_bereich_min[i] = flux_bereich.min()

# Abarbeiten der Spektrenliste:
for i in range(len(filelist)):
    print(i, filelist[i], ' being processed')
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']
    JD_float = float(JD)
    wave_bereich = np.array([])
    flux_bereich = np.array([])

    if (JD_float >= JD_anfang) and (JD_float <= JD_ende):
        wave = np.zeros(header['NAXIS1'])
        for k in range(header['NAXIS1']):
            wave[k] = lambda0 + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])


    plt.scatter(wave_bereich, np.full(len(flux_bereich), JD),
                c=flux_bereich, s=1., cmap='seismic') # adjust the width s of the spectra strips
    plt.clim(flux_bereich_min.min()*0.95, flux_bereich_max.max()*1.05)
    plt.pause(.05)

plt.figure(1)
plt.clim(flux_bereich_min.min()*0.95, flux_bereich_max.max()*1.05)
plt.colorbar()

plt.savefig("LineDevelopment"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.pdf', format='pdf')
plt.savefig("LineDevelopment"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.png', format='png')

print('End of the program')
