#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script is used to model an absorption or emission line profile.
The profile can have one or two maxima or minima.

Various models (Gauss, Voigt, SplitLorentz for asymmetric profiles) can be selected.
can be selected.

The script loads a series of spectra normalized to the continuum in fits format
format. If the line profile in a series is Doppler shifted or completely different from the first
from the first spectrum, the fitting may not work.
Such spectra must then be treated individually.

The first spectrum is plotted and the line to be modeled can be
for a time defined and changeable in line 89 either by clicking with the mouse
and select it or alternatively manually change the wavelength range and the wavelength of the
wavelength of the line extremum.

When graphically selecting the line, you must click four times with the mouse in the plot
first to the left of the line (with as much continuum as possible), then to the right
(with as much continuum as possible) and thirdly the approximate stronger
line minimum/maximum and as fourth the weaker line minimum/maximum
(or again the line minimum/maximum).
This selects the relevant wavelengths and fluxes.

The model to be used for fitting can be selected. The following are available
Gauss, Doppelgauss, Lorentz, Voigt, Doppelvoigt.
The fitting results are saved in an Excel file.

The graphics can be saved if required.

Translated with DeepL.com (free version)

20240326

@author: lothar schanne
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from lmfit.models import LinearModel, GaussianModel, SplitLorentzianModel, VoigtModel
import glob
import pandas as pd


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend
plt.ion()

# Create file list. Spectra in a (sub)folder.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a chronological order
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste:")
print(filelist)
print("\nAnzahl der Spektren: ", len(filelist), "\n")

# Grafik erstes Spektrum
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

# Possibly adjust the waiting time in which you can increase the spectrum:
plt.pause(10)


frage1 = input(
    "Would you like to enter the spectrum limits and the line minima/maxima \
numbers (m) or by mouse click (graphically, g)? Enter m or g: "
)

if frage1 == "g":
    # Interactive definition of the wavelength limits left side with as much
    # continuum as possible, then right side with as much continuum as possible,
    # then click on the stronger minimum/maximum of the line, then click on
    # the weaker second minimum/maximum.
    pts = []
    pts = np.asarray(plt.ginput(n=4, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    extremum1 = pts[2, 0]
    extremum2 = pts[3, 0]
    print("Selected range:", begin, " to", end, ', Extremum1 =', extremum1,
          'Extremum2 =', extremum2)
    print()

if frage1 == "m":
    begin = float(
        input("Enter the short-wave wavelength limit: ")
    )
    end = float(
        input("Enter the long-wave wavelength limit: "))
    print()
    extremum1 = float(
        input("Enter the approximate wavelength of the stronger line minimum/maximum: "))
    extremum2 = float(
        input("Enter the approximate wavelength of the weaker line minimum/maximum: "))

plt.ioff()

modell = input('Which model do you want to use? \
 Gauss (1), double-Gauss (2), Split-Lorentzian (3), Voigt (4), double-Voigt (5)\
 Please select the indicated number: ')

graphspeichern = input('Would you like to save the graphics? Then enter y:')

# Processing the filelist
for k in np.arange(len(filelist)):
    print()
    print(filelist[k])
    sp = fits.open(filelist[k], ignore_missing_end=True)
    # Readout of the observation time as JD, header check
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


# Model definitions:

    # Model for the continuum:
    lin_mod = LinearModel(prefix='lin_')
    pars = lin_mod.make_params(intercept=1, slope=0)

    # Models for the line:
    if modell == '1':
        # Gauss:
        gauss1 = GaussianModel(prefix='g1_')
        if pts[2,1] > 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))

        mod = lin_mod + gauss1

    if modell == '2':
        # Double gauss if the second gauss line overlaps:
        gauss1 = GaussianModel(prefix='g1_')
        if pts[2,1] > 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(gauss1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        gauss2 = GaussianModel(prefix='g2_')
        if pts[3,1] > 1:
           pars.update(gauss2.make_params(center=dict(value=extremum1),
                                      sigma=dict(value=2, min=0),
                                      amplitude=dict(value=pts[3,1], min=0.)))
        if pts[3,1] < 1:
           pars.update(gauss2.make_params(center=dict(value=extremum1),
                                      sigma=dict(value=2, min=0),
                                      amplitude=dict(value=pts[3,1]-1., max=0.)))
        mod = lin_mod + gauss1 + gauss2

    if modell == '3':
        # SplitLorentzian model, designed for asymmetrical lines
        lorentz_mod = SplitLorentzianModel(prefix='lor_')
        if pts[2,1] > 1:
            pars.update(lorentz_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       sigma_r=dict(value=10, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(lorentz_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       sigma_r=dict(value=10, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        mod = lin_mod + lorentz_mod

    if modell == '4':
        # Voigt-Model
        voigt_mod = VoigtModel(prefix='voi_')
        if pts[2,1] > 1:
            pars.update(voigt_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(voigt_mod.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=5, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        mod = lin_mod + voigt_mod

    if modell == '5':
        # double Voigt-Model
        voigt_mod1 = VoigtModel(prefix='voi_1')
        if pts[2,1] > 1:
            pars.update(voigt_mod1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=3, min=0),
                                       amplitude=dict(value=pts[2,1], min=0.)))
        if pts[2,1] < 1:
            pars.update(voigt_mod1.make_params(center=dict(value=extremum1),
                                       sigma=dict(value=3, min=0),
                                       amplitude=dict(value=pts[2,1]-1., max=0.)))
        voigt_mod2 = VoigtModel(prefix='voi_2')
        if pts[3,1] > 1:
            pars.update(voigt_mod2.make_params(center=dict(value=extremum2),
                                       sigma=dict(value=1, min=0),
                                       amplitude=dict(value=pts[3,1], min=0.)))
        if pts[3,1] < 1:
            pars.update(voigt_mod2.make_params(center=dict(value=extremum2),
                                       sigma=dict(value=1, min=0),
                                       amplitude=dict(value=pts[3,1]-1., max=0.)))
        mod = lin_mod + voigt_mod1 + voigt_mod2

    # Fit the model
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

    # Plotten
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    axes[0].plot(x, y)
    axes[0].plot(x, init, '--', label='initial fit')
    axes[0].plot(x, out.best_fit, '-', label='best fit')
    axes[0].legend()
    axes[0].set_xlabel("Wavelenght in Angström", fontsize=10)

    comps = out.eval_components(x=x)
    axes[1].plot(x, y)
    axes[1].plot(x, comps['lin_'], '--', label='Continuum Component')
    if modell == '1':
        axes[1].plot(x, comps['g1_'], '--', label='Gauss Component 1')
    if modell == '2':
        axes[1].plot(x, comps['g1_'], '--', label='Gauss Component 1')
        axes[1].plot(x, comps['g2_'], '--', label='Gauss Component 2')
    if modell == '3':
        axes[1].plot(x, comps['lor_'], '--', label='Lorentz Component')
    if modell == '4':
        axes[1].plot(x, comps['voi_'], '--', label='Voigt Component')
    if modell == '5':
        axes[1].plot(x, comps['voi_1'], '--', label='Voigt Component 1')
        axes[1].plot(x, comps['voi_2'], '--', label='Voigt Component 2')
    axes[1].legend()
    axes[1].set_xlabel("Wavelenght in Angström", fontsize=10)
    plt.pause(.2)

    if graphspeichern == 'y':
        plt.savefig('Linefit_'+filelist[k] + '.pdf')

plt.show()


# Save fitted model parameters as an Excel file in the working directory
ergebnis.to_excel('Fitting_Results.xlsx', sheet_name='Results')


fr = input('When you have finished viewing the graphics and wish to exit the \
program, press the e button: ')
if fr == 'e':
    print('End of ptogram')
    plt.close('all')
