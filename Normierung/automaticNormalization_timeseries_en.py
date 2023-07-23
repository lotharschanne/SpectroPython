#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Works for a single or a time series of 1d spectra in fits format.
At first spectrum the parameters 'inter' and 'sm' are optimized.
Takes from pixel to pixel an interval of 'inter' pixels, calculated
of which is a percentile and compares it with the neighboring intervals.
Finds so local maxima = support points.
You can also distinguish several wavelength ranges from the formation of
normalization points, i.e. delete them (wide lines).
The Smoothing parameter 'sm' adjusts the sensitivity of the spline
which is used for normalization at the end. After that, you can still delete
points that are a certain percentage larger or smaller than the spline.
Finally, all spectra of the time series are normalized with the optimized
parameters. Normalized spectra (and the plots and the normalization functions)
will be saved.

07.12.2020
@author: lothar
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import sigmaclip
import glob
import concurrent.futures as futures
import os
from csaps import csaps

from numba import jit


# definition of functions

# @jit
def wavecalc(flux, header):
    flux_median = np.nanmedian(flux)
    wave = np.zeros(header["NAXIS1"])
    if "CRPIX1" not in header:
        header["CRPIX1"] = 1
    header["CRVAL1"] = header["CRVAL1"] + \
        (1 - header["CRPIX1"]) * header["CDELT1"]
    for i in np.arange(header["NAXIS1"]):
        wave[i] = header["CRVAL1"] + i * header["CDELT1"]
    wave = wave[flux >= flux_median * 0.001]
    flux = flux[flux >= flux_median * 0.001]
    return wave, flux


# @jit
def maximacalc(inter, flux, wave):
    pixelanzahl = len(wave)
    maximas_flux = np.zeros(pixelanzahl)
    maximas_wave = np.zeros(pixelanzahl)
    points_flux = np.zeros(pixelanzahl)
    points_wave = np.zeros(pixelanzahl)
    for n in np.arange(inter, pixelanzahl - 2 * inter):
        interflux = flux[n: n + inter]  # Flux des Intervalls n
        inter_percentil = np.percentile(interflux, 50)
        b = np.percentile(flux[n - inter: n], 75)
        c = np.percentile(flux[n + inter: n + 2 * inter], 75)
        if b < inter_percentil > c:
            interflux1, low, up = sigmaclip(interflux, 2, 2)
            maximas_flux[n + inter // 2] = np.percentile(interflux1, 50)
            maximas_wave[n + inter // 2] = wave[n + inter // 2]
    points_flux = maximas_flux[maximas_flux != 0]
    points_wave = maximas_wave[maximas_wave != 0]
    return points_flux, points_wave


@jit
def reduction(s1, s2, points_wave, points_flux, wave, normfunction):
    for n in np.arange(len(points_wave)):
        for i in np.arange(len(wave)):
            if wave[i] == points_wave[n] and (
                (
                    points_flux[n] < (1 - s1 / 100) * normfunction[i]
                    or points_flux[n] > (1 + s2 / 100) * normfunction[i]
                )
            ):
                points_flux[n] = 0
                break
    return points_flux


@jit
def cutting(wave, points_wave):
    for n in np.arange(len(wave)):
        if wave[n] < points_wave[0]:
            wave[n] = 0
        if wave[n] > points_wave[-1]:
            wave[n] = 0
    return wave


# @jit
def pipe(kl):
    flux, header = fits.getdata(filelist[kl], header=True)
    # Generation of arrays with the wavelengths and fluxes of the spectrum
    wave, flux = wavecalc(flux, header)

    points_flux, points_wave = maximacalc(inter, flux, wave)

    for i in np.arange(len(loesch1)):
        for n in np.arange(len(points_wave)):
            if loesch1[i] < points_wave[n] and loesch2[i] > points_wave[n]:
                points_flux[n] = 0

    points_wave = points_wave[points_flux != 0]
    points_flux = points_flux[points_flux != 0]

    # Restriction of the wavelength range to the range of the interpolation points
    wave = cutting(wave, points_wave)

    flux = flux[wave != 0]
    wave = wave[wave != 0]

    # Creation of a cubic spline from the points (= normfunction)
    normfunction = csaps(points_wave[:], points_flux[:], wave[:], smooth=sm)

    if frage2[0] == "y":
        for m in range(len(s1)):
            points_flux = reduction(
                s1[m], s2[m], points_wave, points_flux, wave, normfunction
            )
            points_wave = points_wave[points_flux != 0]
            points_flux = points_flux[points_flux != 0]
            wave = cutting(wave, points_wave)
            flux = flux[wave != 0]
            wave = wave[wave != 0]
            normfunction = csaps(
                points_wave[:], points_flux[:], wave[:], smooth=sm)

    # Normalization
    normspectrumflux = flux / normfunction
    ausdruck = (
        "Spectrum "
        + str(kl + 1)
        + " "
        + filelist[kl]
        + " Number of points "
        + str(len(points_wave))
    )
    print(ausdruck)

    return (
        kl,
        normfunction,
        normspectrumflux,
        points_flux,
        points_wave,
        wave,
        flux,
        header,
    )


######
# Initial values of changeable constants, later adapted to the line density
# and line width and number of interpolation points of spline
inter = 10  # Intervall (Pixel)

sm = 0.001  # initial smoothing factor for spline

#####

# plt.ion()

# ++++++++++++++
# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results # in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of spectra:")
print(filelist)
print("\nNumber of spectra: ", len(filelist), "\n")


# ********** Optimizing the parameters on the first spectrum
print(
    'We first optimize the parameter  "interval width" and select wavelength ' +
    'ranges in which no interpolation points are desired. Please wait. \n'
)

k = 0
flux, header = fits.getdata(filelist[0], header=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

# Generation of arrays with the wavelengths and fluxes of the spectrum
wave, flux = wavecalc(flux, header)
# The list wave contains the wavelengths of the pixels.
# In the list flux the corresponding intensities.

# Plot the spectrum
fig, ax = plt.subplots(1, figsize=(10, 8))
ax.plot(wave, flux)
ax.set_xlabel("Wavelength [Angström]")
ax.set_ylabel("ADU")
fig.suptitle("Spectrum " + filelist[0])
ax.grid(True)


points_flux, points_wave = maximacalc(inter, flux, wave)

# Plot of normalization points
ax.plot(points_wave, points_flux, "or", markersize=3)
plt.show(block=False)

# +++++++++++++++++++ Intervall width adjusting
print('Current interval width "inter" is:', inter)
print("\nCurrent number of support points: ", len(points_wave))
forderung = input(
    'Do you want to change the interval width "inter"? Then enter y: ')
while forderung == "y":
    inter = int(input("Enter the new interval width: "))
    fig = plt.figure(figsize=(10, 8))
    plt.plot(wave, flux)
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("ADU")
    plt.title("Spectrum " + filelist[0])
    plt.grid(True)

    points_flux, points_wave = maximacalc(inter, flux, wave)

    # Plot of normalization points
    plt.plot(points_wave, points_flux, "or", markersize=3)
    plt.show(block=False)

    print("Current number of support points: ", len(points_wave))

    forderung = input(
        'Do you want to change the interval width "inter" again? Then enter y: '
    )


# ++++++++++++++++++ Exclude wavelength ranges
forderung = input(
    "\nDo you want to remove interpolation points in certain wavelength ranges? Then enter y: "
)
loesch1 = np.zeros(100, dtype=float)
loesch2 = np.zeros(100, dtype=float)

i = 0

while forderung == "y":
    loesch1[i] = float(input("Initial wavelength of the range? :"))
    loesch2[i] = float(input("End wavelength of the range? :"))
    for n in np.arange(len(points_wave)):
        if loesch1[i] < points_wave[n] and loesch2[i] > points_wave[n]:
            points_flux[n] = 0
    i += 1
    forderung = input(
        "Do you want to remove further support points? Then enter y: ")

loesch1 = loesch1[loesch1 != 0]
loesch2 = loesch2[loesch2 != 0]

points_wave = points_wave[points_flux != 0]
points_flux = points_flux[points_flux != 0]

wave = cutting(wave, points_wave)

flux = flux[wave != 0]
wave = wave[wave != 0]
print("Current number of support points: ", len(points_wave))

fig = plt.figure(figsize=(10, 8))
plt.xlabel("Wavelength [Angström]")
plt.ylabel("normalized flux")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
# plt.ylim(flux.min(), flux.max())
plt.plot(wave, flux)
plt.plot(points_wave, points_flux, "or", markersize=3)


# Creation of a cubic spline from the points (= normfunction)
normfunction = csaps(points_wave[:], points_flux[:], wave[:], smooth=sm)
plt.plot(wave, normfunction, "k")
plt.show(block=False)


# +++++++++++++++++++++ Change Smoothing Factor
print('\nCurrently the smoothing factor is "sm": ', sm)
print(
    'With the Smoothing factor we influence the "detail" of the spline: \nsm = 1 means that the spline goes through all points, \nsm = 0 means linear regression through all points (straight line). \nSelect a number between 0. and 1.'
)

sm_Frage = input(
    "\nDo you want to change the smoothing factor? Then enter y: ")

while sm_Frage == "y":
    sm = float(input("Enter the new value for sm: "))
    normfunction = csaps(points_wave[:], points_flux[:], wave[:], smooth=sm)
    fig = plt.figure(figsize=(10, 8))
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("normalized flux")
    plt.title("Spectrum " + filelist[0])
    plt.grid(True)
    plt.plot(wave, flux)
    plt.plot(points_wave, points_flux, "or", markersize=3)
    plt.plot(wave, normfunction, "k")
    plt.show(block=False)
    sm_Frage = input(
        "Do you want to change the smoothing factor again? Then enter y: ")


# ************** Remove interpolation points below and above the spline
s1 = []
s2 = []
frage2 = []
j = 0
frage2.append(
    input(
        "\nIf you want to delete interpolation points s1 % below and s2 % above the spline, enter y: "
    )
)
while frage2[j] == "y":
    s1.append(float(input("Enter s1 in %: ")))
    s2.append(float(input("Enter s2 in %: ")))
    print("Please wait")
    points_flux = reduction(s1[j], s2[j], points_wave,
                            points_flux, wave, normfunction)

    points_wave = points_wave[points_flux != 0]
    points_flux = points_flux[points_flux != 0]

    wave = cutting(wave, points_wave)

    flux = flux[wave != 0]
    wave = wave[wave != 0]
    print("Current number of support points: ", len(points_wave))

    normfunction = csaps(points_wave[:], points_flux[:], wave[:], smooth=sm)

    fig = plt.figure(figsize=(10, 8))
    plt.xlabel("Wavelength [Angström]")
    plt.ylabel("normalized flux")
    plt.title("Spectrum " + filelist[0])
    plt.grid(True)
    plt.plot(wave, flux)
    plt.plot(points_wave, points_flux, "or", markersize=3)
    plt.plot(wave, normfunction, "k")
    plt.show(block=False)

    j += 1
    frage2.append(
        input(
            "\nIf you want to delete further interpolation points s1 % below and s2 % above the spline, enter y: "
        )
    )


normspectrumflux = flux / normfunction

fig = plt.figure(figsize=(10, 8))
plt.xlabel("Wavelength [Angström]")
plt.ylabel("normalized flux")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
plt.ylim(flux.min(), flux.max())
plt.plot(wave, flux)
plt.plot(points_wave, points_flux, "or", markersize=3)
plt.plot(wave, normfunction, "k")
plt.show(block=False)

frage3 = input("If the graphic is to be saved as pdf, then enter y: ")
if frage3 == "y":
    filename = filelist[0].rsplit(".", 1)
    filename = filename[0] + "_points.pdf"
    plt.savefig(filename)


# Plot the normalized spectrum
fig = plt.figure(figsize=(10, 8))
plt.plot(wave, normspectrumflux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("normalized flux")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
plt.show(block=False)

if frage3 == "y":
    filename = filelist[0].rsplit(".", 1)
    filename = filename[0] + "_normalized.pdf"
    plt.savefig(filename)


# saving the normalized spectrum as fit
antwort = input(
    "\nIf you want to save the normalized spectrum as fit, then enter y: ")
if antwort == "y":
    header["CRVAL1"] = wave[0]
    header["NAXIS1"] = len(wave)
    header["CRPIX1"] = 1
    filename = filelist[0].rsplit(".", 1)
    filename = filename[0] + "_normalized.fit"
    fits.writeto(
        filename, normspectrumflux, header, overwrite=True, output_verify="silentfix"
    )

# saving the normfunction
antwort1 = input("If you want to save the normfunction as fit, then enter y: ")
if antwort1 == "y":
    header["CRVAL1"] = wave[0]
    header["NAXIS1"] = len(wave)
    header["CRPIX1"] = 1
    filename = filelist[0].rsplit(".", 1)
    filename = filename[0] + "_normfunction.fit"
    fits.writeto(filename, normfunction, header, overwrite=True)


# *********** Normalize entire time series  *************

frage1 = input(
    "\nIf you want to normalize the other spectra of the time series, enter y: "
)
print(
    "If your time series contains many (> 10) spectra, the displayed graphs of raw and normalized spectra may block a lot of memory."
)
frage4 = input(
    "Do you want to calculate and to see all plots with the raw spectra, normalization points and the normalization function? Then enter y: "
)
frage5 = input(
    "Do you want to calculate and to see all plots with the normalized spectra? Then enter y: "
)
print(
    "\nNow please wait until the program is finished. This may take some time. With Ctrl C it can be aborted. The normalized spectra and pdf of the diagrams can be found in the working folder."
)

if __name__ == "__main__":
    if frage1 == "y":
        with futures.ProcessPoolExecutor(max_workers=os.cpu_count() - 1) as e:
            fs = {e.submit(pipe, n): n for n in np.arange(1, len(filelist))}
            fertig = False
            while not fertig:
                res = futures.wait(fs)
                for f in res.done:
                    (
                        k,
                        normfunction,
                        normspectrumflux,
                        points_flux,
                        points_wave,
                        wave,
                        flux,
                        header,
                    ) = f.result()
                    del fs[f]
                    # Saving the normalized spectrum as fit
                    if antwort == "y":
                        filename = filelist[k].rsplit(".", 1)
                        filename = filename[0] + "_normalized.fit"
                        header["CRVAL1"] = wave[0]
                        header["NAXIS1"] = len(wave)
                        header["CRPIX1"] = 1
                        fits.writeto(
                            filename,
                            normspectrumflux,
                            header,
                            overwrite=True,
                            output_verify="silentfix",
                        )
                    # saving the normfunction
                    if antwort1 == "y":
                        header["CRVAL1"] = wave[0]
                        header["NAXIS1"] = len(wave)
                        header["CRPIX1"] = 1
                        filename = filelist[k].rsplit(".", 1)
                        filename = filename[0] + "_normfunction.fit"
                        fits.writeto(filename, normfunction,
                                     header, overwrite=True)
                    # plotting
                    if frage4 == "y":
                        # plot
                        fig = plt.figure(figsize=(10, 8))
                        plt.ylim(flux.min(), flux.max())
                        plt.plot(wave, flux)
                        plt.xlabel("Wavelength [Angström]")
                        plt.ylabel("ADU")
                        plt.title("Spectrum " + filelist[k])
                        plt.grid(True)
                        plt.plot(points_wave, points_flux, "or", markersize=3)
                        plt.plot(wave, normfunction, "k")
                        if frage3 == "y":
                            filename = filelist[k].rsplit(".", 1)
                            filename = filename[0] + "_points.pdf"
                            plt.savefig(filename)
                        plt.show(block=False)

                    if frage5 == "y":
                        # Plot the normalized spectrum
                        fig = plt.figure(figsize=(10, 8))
                        plt.plot(wave, normspectrumflux)
                        plt.xlabel("Wavelength [Angström]")
                        plt.ylabel("normalized flux")
                        plt.title("Spectrum " + filelist[k])
                        plt.grid(True)
                        if frage3 == "y":
                            filename = filelist[k].rsplit(".", 1)
                            filename = filename[0] + "_normalized.pdf"
                            plt.savefig(filename)
                        plt.show(block=False)

                fertig = len(res.not_done) == 0

x = input("Type e when you are finished viewing the graphics: ")
while x:
    if x == "e":
        print("END of program, all done")
        x = False
    else:
        x = input(
            "The input was not e.Type e when you are finished viewing the graphics: "
        )
