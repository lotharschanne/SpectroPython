#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 22:01:37 2023

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Eingabe des Bereichs der Beobachtungszeitpunkte als JD
JD_anfang = float(input('Geben Sie den Anfang des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))
JD_ende = float(input('Geben Sie das Ende des Bereichs der \
Beobachtungszeitpunkte (JD) ein: '))

# Eingabe des abgebildeten Wellenlängenbereichs
lambda_anfang = float(input('Geben Sie den Anfang des abzubildenden \
Wellenlängenbereichs an: '))
lambda_ende = float(input('Geben Sie das Ende des abzubildenden \
Wellenlängenbereichs an: '))

obj = input("Please enter the object name: ")

plt.figure(1, figsize=(7, 10))
# Abarbeiten der Spektrenliste:
for i in range(len(filelist)):
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

    if JD_float >= JD_anfang and JD_float <= JD_ende:
        wave = np.zeros(header['NAXIS1'])
        for k in range(header['NAXIS1']):
            wave[k] = lambda0 + k * step
            if wave[k] >= lambda_anfang and wave[k] <= lambda_ende:
                wave_bereich = np.hstack([wave_bereich, wave[k]])
                flux_bereich = np.hstack([flux_bereich, flux[k]])

    plt.scatter(np.full(len(flux_bereich), JD), wave_bereich,
                c=flux_bereich, s=0.5, cmap='cubehelix')

    plt.clim(1., 2.)  # anpassen an auftrende Fluxe
    plt.pause(.05)

plt.colorbar()

plt.title("Linienentwicklung" + obj)
plt.grid(True)
plt.xlabel('Wellenlänge [Angström]')
plt.savefig("Linienentwicklung"+'JD_'+str(JD_anfang)+'_'+str(JD_ende)
            + '_'+str(lambda_anfang)+'_'+str(lambda_ende) + '.pdf', format='pdf')
