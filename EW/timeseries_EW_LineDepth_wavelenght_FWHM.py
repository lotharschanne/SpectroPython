#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculates for a time series of ON THE CONTINUUM NORMALIZED SPECTRES in the
fits format, the equivalent width EW of a line, the line depth, the FWHM and
wavelength of the line minimum.
Enter the integration limits manually or graphically.
Normalization errors are compensated by a linear renormalization routine in the
integration range. The flux at the integration limits is taken for the continuum
calculation and linearly interpolated between this 2 points
and normalized again with the resulting linear continuum function.
It is therefore important that the integration limits actually lie on the continuum.

The calculations assume that the same wavelength interval is used for all
spectra in the series and no significant RV changes take place, or if they do,
that the line is isolated (i.e. the flux = 1 in the environment). Therefore best
use barycentrically corrected spectra.

The first spectrum is plotted and displayed for 15 seconds. During this
period, the graphics window is interactive, so that you can enlarge the spectrum
and the integration wavelength range for the EW calculation can be selected optically.
This reaction time of 20 seconds can be changed in line 92. If the graphical/
interactive input of the integration limits was selected, the integration range
must be selected after the pause with two mouse clicks in the graphics window.
Otherwise the integration limits must be entered as numbers.

The selected line is plotted as a black line for all spectra,
the renormalized line in blue and the interpolated line in red.
On request, the plots can also be saved as a pdf file.

Version 2023-05-08

@author: lothar schanne
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import Rbf

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Question whether the graphics should be saved
frage = input('Should the graphics of the spectra sections be saved as a PDF? \
If yes, please enter "j", otherwise enter another button or return: ')

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout the list
print("\nSpectrum list: \n")
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

frage1 = input(
    "Would you like to enter the integration limits numerically or \
by mouse click (graphically)? Enter m or g: "
)

# Plot the spectrum
fig = plt.figure()
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("ADU")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
print('\nPlease narrow the plot interactively to the integration area within 20 seconds')
plt.pause(15)
print('Now set the integration range per 2 mouse clicks')


if frage1 == "g":
    # Interactive definition of the integration limits
    # If the display of the spectrum is interactively enlarged, delete the input points.
    # Delete them with the right mouse button until the left
    # side of the line to be measured is clicked.
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    print("Selected integration range:", begin, " bis", end)
    print()

if frage1 == "m":
    # Eingabe der Integrationsgrenzen
    begin = float(
        input("Enter the short-wave wavelength integration limit: ")
    )
    end = float(
        input("Enter the long-wave wavelength integration limit: "))
    print()

# Initiating arrays
EW = np.zeros(len(filelist))
JD = np.zeros(len(filelist))
linienminimum_flux = np.zeros(len(filelist))
linienminimum_wave = np.zeros(len(filelist))
fwhm = np.zeros(len(filelist))

# Processing the filelist
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)

    # Header Überprüfung
    if "JD-OBS" in header:
        JD[i] = float(header["JD-OBS"])
    elif 'JD' in header:
        JD[i] = float(header["JD"])
    elif "JD_OBS" in header:
        JD[i] = float(header["JD_OBS"])
    elif "BAS_MJD" in header:
        JD[i] = float(header["BAS_MJD"])
    else:
        print("No observation time in the header")
        break

    try:
        header["CRPIX1"]
    except:
        header["CRPIX1"] = 1

    # Generation of array with the wavelengths of the spectrum
    wave = np.ones(len(flux), dtype=float)
    for k in range(header["NAXIS1"]):
        wave[k] = (
            header["CRVAL1"]
            + (k - header["CRPIX1"] + 1) * header["CDELT1"])

    # Isolate the integration range:
    ind = []
    for n in range(len(wave)):
        if wave[n] >= begin and wave[n] <= end:
            ind.append(n)
    intervallwave = wave[ind]
    intervallflux = flux[ind]

    fig = plt.figure()
    plt.plot(intervallwave, intervallflux, '-k',
             label='Original measurement')

    intervallflux_begin = intervallflux[0]
    intervallflux_end = intervallflux[-1]
    # Slope calculation for linear interpolation of the continuum
    faktor = (intervallflux_end - intervallflux_begin) / len(intervallflux)

    # Renormalization of the integration domain to the (linear) pseudo-continuum
    # formed from begin and end (the integration limit wavelengths)
    ew = 0
    for p in range(len(intervallflux)):
        intervallflux[p] = intervallflux[p] / \
            (intervallflux_begin + p * faktor)
        # Berechnung der EW
        dif = 1 - intervallflux[p]
        ew += dif * header["CDELT1"]
    EW[i] = ew

    # Increasing the resolution tenfold through interpolation
    # Radial basis function (RBF) over the spectrum in the interval
    # Smooth adjustment, 0. = function goes through all points, >0. = equalization
    rbf = Rbf(intervallwave, intervallflux, smooth=1)  # adjust smooth
    newwaveintervall = np.arange(
        intervallwave[0], intervallwave[-1], header['CDELT1'] / 10)
    newfluxinterpolated = rbf(newwaveintervall)
    linienminimum_flux[i] = newfluxinterpolated.min()
    linienminimum_wave[i] = newwaveintervall[newfluxinterpolated.argmin()]
    print(filelist[i] + ':')
    print('EW:', EW[i])
    print('Wavelength Line minimum = {:.2f}'.format(
        linienminimum_wave[i]), 'Angström')
    print('Flux line minimum = {:.3f}'.format(linienminimum_flux[i]))

    # FWHM ermitteln
    hm = (1 + linienminimum_flux[i]) / 2
    for o in range(newfluxinterpolated.argmin(), 0, -1):
        if newfluxinterpolated[o] >= hm:
            links = newwaveintervall[o]
            break
    for o in range(newfluxinterpolated.argmin(), len(newfluxinterpolated), 1):
        if newfluxinterpolated[o] >= hm:
            rechts = newwaveintervall[o]
            break
    fwhm[i] = rechts - links
    print('FWHM = {:.2f}'.format(fwhm[i]), ' Angström\n')

    # Plotten der Linien
    plt.plot(intervallwave, intervallflux, color='b',
             label='renormalized spectrum',
             linewidth=1.)
    plt.plot(linienminimum_wave[i],
             linienminimum_flux[i], "or",
             label='Minimum')
    plt.title('Used spectrum range from ' + filelist[i])
    plt.plot(newwaveintervall, newfluxinterpolated, color='r', linewidth=.6,
             label='interpolated spectrum')
    plt.legend(loc='best')
    plt.grid(True)
    plt.pause(.1)

    if frage == 'j':
        plt.savefig(filelist[i] + '{:.2f}_'.format(begin) + '{:.2f}'.format(end)
                    + '.pdf')

    plt.close('all')

# Saving the calculated data
ascii.write(
    [filelist, JD, EW, linienminimum_wave, linienminimum_flux, fwhm],
    'Wavelength_range' +
    '{:.2f}_'.format(begin) + '{:.2f}'.format(end) + '.csv',
    names=["Spectrum ", 'JD', "EW", 'Wave length line minimum',
           'Flux Line minimum', 'FWHM'],
    overwrite=True,
    format="csv",
)


# # Plotting
# fig = plt.figure()
# plt.stem(JD, -EW)
# plt.xlabel("JD")
# plt.ylabel("-EW in Angstroem")
# plt.title('EW ' + str(begin) + ' bis ' + str(end) + ' Angstroem')
# plt.grid(True)
# plt.savefig('EW_' + str(begin) + '_' + str(end) + '.pdf', format='pdf)
# plt.savefig('EW_' + str(begin) + '_' + str(end) + '.png', format='png)

print('Ende des Programms')
