#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Das Skript liest eine Serie von fits-1d-Spektren ein und berechnet aus
unterschiedlichen Headereinträgen für das Beobachtungsdatum das jeweilige JD
und trägt es in den header des jeweiligen Spektrums ein.

Created on Thu Dec  8 13:35:37 2022

@author: lothar
"""

import glob
from astropy.io import fits
from astropy.time import Time


# filelist erstellen
# Pfad und Name der Spektren eingeben
files = input("Geben Sie Pfad und Namen zu den Spektren ein: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Abarbeitung der Spektrenliste
for i in range(len(filelist)):
    g = fits.open(filelist[i])
    print()
    print(filelist[i])
    if 'JD' not in g[0].header:
        if 'OBSDATE' in g[0].header:
            try:
                t = Time(g[0].header['OBSDATE'])
                g[0].header['JD'] = t.jd
                print('OBSDATE gefunden und als JD in den Header eingetragen')
            except:
                print(filelist[i],
                      ': OBSDATE in falschem Format')
        elif 'JD-OBS' in g[0].header:
            try:
                t = g[0].header['JD-OBS']
                g[0].header['JD'] = t
                print('JD-OBS gefunden und als JD in den Header eingetragen')
            except:
                print(filelist[i],
                      ': JD-OBS in falschem Format')
        elif 'JD-MID' in g[0].header:
            try:
                t = g[0].header['JD-MID']
                g[0].header['JD'] = t + 2400000.5
                print('JD-MID gefunden und als JD in den Header eingetragen')
            except:
                print(filelist[i],
                      ': JD-MID in falschem Format')
        elif 'DATE-OBS' in g[0].header:
            try:
                t = Time(g[0].header['DATE-OBS'])
                g[0].header['JD'] = t.jd
                print('DATE-OBS gefunden und als JD in den Header eingetragen')
            except:
                print(
                    filelist[i], ': DATE-OBS in falschem Format')
        elif 'DATE-BEG' in g[0].header:
            try:
                t = Time(g[0].header['DATE-BEG'])
                g[0].header['JD'] = t.jd
                print('DATE-BEG gefunden und als JD in den Header eingetragen')
            except:
                print(
                    filelist[i], ': DATE-BEG in falschem Format')
        elif "MJD-OBS" in g[0].header:
            try:
                t = g[0].header["MJD-OBS"]
                g[0].header['JD'] = t + 2400000.5
                print('MJD-OBS gefunden und als JD in den Header eingetragen')
            except:
                print(
                    filelist[i], ': MJD-OBS in falschem Format')
        elif "BAS-MJD" in g[0].header:
            try:
                t = g[0].header["BAS-MJD"]
                g[0].header['JD'] = t + 2400000.5
                print('BAS-MJD gefunden und als JD in den Header eingetragen')
            except:
                print(
                    filelist[i], ': BAS-MJD in falschem Format')
        fits.update(filelist[i], g[0].data, g[0].header)
    else:
        print(filelist[i], ' Headereintrag JD bereits vorhanden')
    g.close()
