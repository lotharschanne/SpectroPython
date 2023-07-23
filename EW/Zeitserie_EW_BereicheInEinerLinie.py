#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine auf das Kontinuum normierter Zeitserie im fits-Format
die Äquivalentweite einer in Bereiche aufgeteilten Linie.
Eingabe der Integrationsgrenzen (in Angström) über eine Liste, die anzupassen ist.
Die EW-Berechnung setzt voraus, dass für alle Spektren der Serie das gleiche
Wellenlängenintervall für die Berechnung des Integrals verwendet werden kann,
also keine wesentlichen RV-Änderungen stattfinden. Deshalb am besten baryzentrisch
korrigierte Spektren benutzen.
Die ermittelten EW's und die Flux-Minima und -Maxima der Linien-Bereiche
werden in einer asci-Datei gespeichert.

Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import pandas as pd

# Liste der Wellenlängen, die als Bereichsgrenzen für die EW-Integration dienen
WL = [6559., 6561., 6563., 6565.]  # bitte entsprechend ändern

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
# print("\nSpektrenliste: \n")
# print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")


EW = np.zeros((len(filelist), len(WL)-1))
Minima = np.full((len(filelist), len(WL)-1), np.inf)
Maxima = np.zeros((len(filelist), len(WL)-1))
JD = np.zeros(len(filelist))
Data = []

# Abarbeiten der filelist
for k in np.arange(len(filelist)):
    sp = fits.open(filelist[k], ignore_missing_end=True)
    try:
        JD[k] = sp[0].header['JD']
    except:
        JD[k] = None
    Data.append((filelist[k], JD[k]))

    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )
    # The list wave contains the wavelengths of the pixels.
    # In the list flux the corresponding intensities.
    # Close the fits-file:
    sp.close()

# Abarbeiten der Integrationsbereiche:
    for m in np.arange(len(WL)-1):
        for n in np.arange(len(wave)):
            if wave[n] >= WL[m] and wave[n] < WL[m+1]:
                EW[k, m] += (1 - flux[n])*sp[0].header['CDELT1']
                if flux[n] < Minima[k, m]:
                    Minima[k, m] = flux[n]
                if flux[n] > Maxima[k, m]:
                    Maxima[k, m] = flux[n]
    for l in np.arange(len(WL)-1):
        Data.append((WL[l], WL[l+1], EW[k, l], Minima[k, l], Maxima[k, l]))


serie = pd.DataFrame(
    Data, columns=['Wave-left', 'Wave-right', 'EW', 'Flux-Minimum', 'Flux-Maximum'])
fileobj = open('EWs.csv', 'w')
fileobj.writelines(serie.to_csv())
fileobj.close()
