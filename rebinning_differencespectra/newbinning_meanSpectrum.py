#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a spectrum catalog (fits), rebinds the spectra and generates tab
spectra (tab table with the column names WAVE and FLUX) and fits files, all
with the same selectable step size and the same wavelength range. The type of
interpolation (linear or by cubic spline) can be selected by commenting it out.
In addition, an averaged spectrum is calculated and saved as fits and as a tab
table with the column names WAVE and FLUX. The averaged spectrum is also plotted.

@author: lothar schanne
24.01.2022
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from specutils import Spectrum1D
from specutils.manipulation import (
    LinearInterpolatedResampler,
    SplineInterpolatedResampler,
)
from astropy import units as u
from astropy.visualization import quantity_support
quantity_support()

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print(filelist, end='\n')
print("\nNumber of spectra: ", len(filelist), "\n")

# Calculation and output of the common wavelength range and step sizes
well_min = 0
well_max = 10000
step_min = 100
step_max = 0
for i in range(len(filelist)):
    sp = fits.open(filelist[i])
    crval = sp[0].header["CRVAL1"]
    if "CRPIX1" not in sp[0].header:
        sp[0].header["CRPIX1"] = 1
    crpix = sp[0].header["CRPIX1"]
    cdel = sp[0].header["CDELT1"]
    wave_erstesPixel = crval - cdel * (crpix - 1)
    if wave_erstesPixel > well_min:
        well_min = wave_erstesPixel
    if wave_erstesPixel + cdel * sp[0].header["NAXIS1"] < well_max:
        well_max = wave_erstesPixel + cdel * sp[0].header["NAXIS1"]
    if cdel < step_min:
        step_min = cdel
    if cdel > step_max:
        step_max = cdel
    sp.close()

print("\nCommon wavelength range: ", well_min, well_max)
print("Minimum and maximum step width: ", step_min, step_max)

# Step size and wavelength range to be selected
newstep = float(input("Enter the desired step width:"))
print("\nSpecification of the wavelength range to be transferred ")
a = float(input("Wavelength Begin: "))
b = float(input("Wavelength End: "))

pixzahl = int((b - a) / newstep)
DATA = np.zeros((len(filelist), pixzahl+1))

# Calculation of the new spectra of the series
for i in range(len(filelist)):
    sp = fits.open(filelist[i])
    crval = sp[0].header["CRVAL1"]
    if "CRPIX1" not in sp[0].header:
        sp[0].header["CRPIX1"] = 1
    crpix = sp[0].header["CRPIX1"]
    cdel = sp[0].header["CDELT1"]
    wave_erstesPixel = crval - cdel * (crpix - 1)
    aindex = int((a - wave_erstesPixel) / cdel)
    bindex = int((b - wave_erstesPixel) / cdel)
    newflux = sp[0].data[aindex:bindex]
    newflux = newflux * u.dimensionless_unscaled
    newwave = np.zeros(bindex - aindex)
    for k in range(len(newwave)):
        newwave[k] = wave_erstesPixel + (aindex + k) * cdel
    newwave = newwave * u.AA

    # new binning:
    input_spec = Spectrum1D(spectral_axis=newwave, flux=newflux)
    new_disp_grid = np.arange(a, b, newstep) * u.AA

    # Interpolation
    # lineare Interpolation:
    linear = LinearInterpolatedResampler()
    new_spec = linear(input_spec, new_disp_grid)
    # Interpolation per Spline:
    # spline = SplineInterpolatedResampler()
    # new_spec = spline(input_spec, new_disp_grid)

    filename = (
        filelist[i].rsplit(".")[0] + "_" + str(int(a)) +
        "_" + str(int(b)) + ".dat"
    )
    # save
    ascii.write(
        [new_spec.spectral_axis, new_spec.flux],
        filename,
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )
    # Calculate and save fits file
    sp[0].header["CRVAL1"] = float(new_spec.spectral_axis[0].value)
    sp[0].header["CRPIX1"] = 1
    sp[0].header["NAXIS1"] = len(new_spec.spectral_axis)
    sp[0].header["CDELT1"] = newstep
    newfile = (
        filelist[i].rsplit(".")[0] + "_" + str(int(a)) +
        "_" + str(int(b)) + ".fit"
    )
    fits.writeto(
        newfile,
        new_spec.flux.value,
        sp[0].header,
        overwrite=True,
        output_verify="silentfix",
    )

    for j in range(len(new_spec.flux)):
        DATA[i, j] = new_spec.flux[j].value
    sp.close()

# Calculating and saving the average spectrum
obj = input("Enter the object name without blanks: ")
mean = np.zeros(len(DATA[0]) - 1)
for m in range(len(mean)):
    mean[m] = DATA[:, m].mean()
file = obj + "_" + str(int(a)) + "_" + str(int(b))
filename = file + "_mean.dat"
ascii.write(
    [new_spec.spectral_axis / u.AA, mean[:]],
    filename,
    overwrite=True,
    names=["WAVE", "FLUX"],
    format="tab",
)
filename = file + "_mean.fits"
fits.writeto(filename, mean[:], sp[0].header,
             overwrite=True, output_verify="silentfix")

# Plotting the mean spectrum
fig = plt.figure(figsize=(14, 14))
plt.xlabel("Wavelength [Angstroem]")
plt.ylabel("Flux")
plt.title("Mittleres Spektrum " + filename)
plt.plot(new_spec.spectral_axis, mean[:])
plt.pause(1)
fig.savefig(filename.rstrip(".")[0] + "_meanSüectrum.png", format='png')
# fig.savefig(filename.rstrip(".")[0] + "_meanSüectrum.pdf", format='pdf')


print('To exit the program, click on the last opened diagram.')
plt.waitforbuttonpress(-1)
plt.close('all')
