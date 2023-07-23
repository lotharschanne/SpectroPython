#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Differenzspektren.py

Liest tab-Spektren einer Serie ein und bildet Differenzspektren zum
angegebenen mittleren Spektrum, die dann als tab-Datei ausgegeben werden.
Alle Differenzspektren werden in einem Plot grafisch dargestellt, der auch
abgespeichert wird.

Created on Sat Feb 27 23:15:12 2021

@author: lothar
"""
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

mean_name = input("Name Mittleres Spektrum: ")

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")
print("Mittleres Spektrum", mean_name)


mean_spec = ascii.read(mean_name, format="tab")

fig = plt.figure(figsize=(14, 14))
plt.xlabel("Wavelength [Angstroem]")
plt.ylabel("Differenz Flux")
plt.grid(True)

for i in range(len(filelist)):
    spec = ascii.read(filelist[i], format="tab")
    diff = np.zeros(len(spec["WAVE"]))
    for j in range(len(spec["WAVE"])):
        diff[j] = spec["FLUX"][j] - mean_spec["FLUX"][j]
    filename = filelist[i].rsplit(".")[0] + "_diffspektrum" + ".dat"
    ascii.write(
        [spec["WAVE"], diff[:]],
        filename,
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )
    plt.plot(spec["WAVE"], diff)

fig.savefig("Differenzspektren" + "overplot.pdf")
plt.show(block=True)
