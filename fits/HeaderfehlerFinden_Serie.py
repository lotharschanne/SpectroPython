#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
HeaderfehlerFinden_Serie.py

Erzeugung einer Spektrenliste der 1d-Spektren im fits-Format in einem Ordner,
Falls im Header ein Fehler ist (Nichtübereinstimmung mit den fits-Kopnventionen),
wird eine Warmeldung ausgegeben.

20230623
@author: lothar
"""

import glob
from astropy.io import fits


# filelist erstellen
# Pfad und Name der Spektren eingeben
files = input("Geben Sie Pfad und Namen zu den Spektren ein: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")


for i in range(len(filelist)):
    g = fits.open(filelist[i])
    g.close(output_verify='warn')
