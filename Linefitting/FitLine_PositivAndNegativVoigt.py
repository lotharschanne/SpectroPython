#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script is used for Voigt modeling of an emission line profile which, as with
beta Lyrae consists of a emission line with a central absorption contained therein.

The script loads a series of spectra normalized to the continuum in fits format
format. If the line profile in a series is Doppler shifted or completely
different from the first spectrum, the fitting may not work.
Such spectra must then be treated individually.

The first spectrum is then plotted and the line to be modeled can be
for a time defined and changeable in line 89 either by mouse click
enlarged and selected.

When graphically selecting the line, you must click four times with the mouse
in the plot, first to the left of the line (with as much continuum as possible),
then to the right (with as much continuum as possible) and thirdly the
approximate line maximum and fourthly on the minimum of the embedded absorption.
This selects the relevant wavelengths and fluxes.

The results of the fitting are saved in an Excel file.

The graphics can be saved if desired.

20240326

@author: lothar schanne
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from lmfit.models import LinearModel, VoigtModel
import glob
import pandas as pd


plt.switch_backend('Qt5Agg')
plt.ion()

# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If named correctly, this results in a chronological order
filelist.sort()

# print of the list
print("\nSpektrenliste:")
print(filelist)
print("\nAnzahl der Spektren: ", len(filelist), "\n")

# Graphic of first spectrum
sp = fits.open(filelist[0], ignore_missing_end=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

if 'CRPIX1' not in sp[0].header:
    sp[0].header['CRPIX1'] = 1

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
print('Now zoom in on the line from the spectrum, then wait')

# Possibly adjust the waiting time in which you can increase the spectrum:
plt.pause(10)

print('Now you can click on the points in the graphic.')

# Interactive definition of the wavelength limits left side with as much
# continuum as possible, then right side with as much continuum as possible,
# then click on the maximum of the line, then click on the absorption minimum.
pts = []
pts = np.asarray(plt.ginput(n=4, timeout=-1))
plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
begin = pts[0, 0]
end = pts[1, 0]
extremum1 = pts[2, 0]
extremum2 = pts[3, 0]
print("Selected range:", begin, " bis", end, ', Maximum =', extremum1,
          'Wavelength of absorption =', extremum2)
print()

plt.ioff()

graphspeichern = input('Would you like to save the graphics? Then enter y:')

# Processing the filelist
for k in np.arange(len(filelist)):
    print()
    print(filelist[k])
    sp = fits.open(filelist[k], ignore_missing_end=True)

    if "JD-OBS" in sp[0].header:
        JD = float(sp[0].header["JD-OBS"])
    elif 'JD' in sp[0].header:
        JD = float(sp[0].header["JD"])
    elif "JD_OBS" in sp[0].header:
        JD = float(sp[0].header["JD_OBS"])
    elif "BAS_MJD" in sp[0].header:
        JD = float(sp[0].header["BAS_MJD"])
    else:
        print("No observation time in the header")
        JD = 0
        break

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

    # Cut out the line area:
    for n in np.arange(len(wave)):
        if wave[n] <= begin and wave[n + 1] > begin:
            begin_n = n
            break
    for n in np.arange(len(wave)):
        if wave[n] <= end and wave[n + 1] > end:
            end_n = n
            break

    x = wave[begin_n:end_n]
    y = flux[begin_n:end_n]


# Definition of models:

    # Model of continuum:
    lin_mod = LinearModel(prefix='lin_')
    pars = lin_mod.make_params(intercept=1., slope=0.)

    # double voigt:
    voigt1 = VoigtModel(prefix='voi_1')
    pars.update(voigt1.make_params(center=dict(value=extremum1),
                                        sigma=dict(value=1., min=0.),
                                        amplitude=dict(value=pts[2,1])))
    voigt2 = VoigtModel(prefix='voi_2')
    pars.update(voigt2.make_params(center=dict(value=extremum2),
                                    sigma=dict(value=3., min=0.),
                                    amplitude=dict(value=-(pts[2,1]-pts[3,1]),
                                                   max=0.)))
    mod = lin_mod + voigt1 + voigt2


    # Fitting of Models
    init = mod.eval(pars, x=x)
    out = mod.fit(y, pars, x=x)

    print(out.fit_report(correl_mode='list'))

    if k == 0:
        ergebnisse = pd.Series(out.values)
        ergebnisse = pd.concat([ergebnisse, pd.Series({'JD': JD})], axis=0)
        ergebnis = pd.DataFrame({filelist[k]: ergebnisse})

    if k > 0:
        ergebnisse = pd.Series(out.values)
        ergebnisse = pd.concat([ergebnisse, pd.Series({'JD': JD})], axis=0)
        ergebnis = pd.concat(
            [ergebnis, pd.DataFrame({filelist[k]:ergebnisse})], axis=1)

    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].plot(x, y)
    axes[0].plot(x, init, '--', label='initial fit')
    axes[0].plot(x, out.best_fit, '-', label='best fit')
    axes[0].legend()
    axes[0].set_xlabel("Wavelength in Angström", fontsize=10)

    comps = out.eval_components(x=x)

    axes[1].plot(x, y)
    axes[1].plot(x, comps['lin_'], '--', label='Continuum Component')
    axes[1].plot(x, comps['voi1_'], '--', label='Voigt Component 1')
    axes[1].plot(x, comps['voi2_'], '--', label='Voigt Component 2')
    axes[1].legend()
    axes[1].set_xlabel("Wavelength in Angström", fontsize=10)
    plt.pause(.2)

    if graphspeichern == 'y':
        plt.savefig('Linefit_'+filelist[k] + '.pdf')

plt.show()


# Gefittete Modellparameter als Excel-Datei im Arbeitsverzeichnis abspeichern
ergebnis.to_excel('Results_Voigt.xlsx', sheet_name='Results')


fr = input('When you have finished viewing the graphics and wish to exit the program, press the e: ')
if fr == 'e':
    print('End of Program')
    plt.close('all')
