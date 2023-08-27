#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")


date_obs = list(range(len(filelist)))

print("Headereinträge")

# Beobachtungszeitpunkt einlesen und ausgeben
print("\nSpektrum", "DATE-OBS")

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True, ignore_missing_end=True)
    if "DATE-OBS" in header:
        date_obs[i] = header["DATE-OBS"]
    elif "DATE_OBS" in header:
        date_obs[i] = header["DATE_OBS"]
    elif "MJD" in header:
        date_obs[i] = header["MJD"]
    elif "JD-MID" in header:
        date_obs[i] = header["JD-MID"]
    else:
    	print("Kein Beobachtunsgzeitpunkt im Header")

    print(filelist[i], date_obs[i])

# Spektraler Bereich berechnen und ausgeben
print("    Spektrumdatei    ", "Anfang   ", "Ende")
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True, ignore_missing_end=True)

    Anfang = header["CRVAL1"]
    Ende = header["CRVAL1"] + header["NAXIS1"] * header["CDELT1"]
    print(filelist[i], f"{Anfang:.1f}, {Ende:.1f}")

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, date_obs],
    "_Beobachtungszeitpunkte" + ".dat",
    overwrite=True,
    names=["Spektrum", "DATE-OBS"],
    format="tab",
)
