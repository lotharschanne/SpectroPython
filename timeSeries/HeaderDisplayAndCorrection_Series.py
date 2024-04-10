#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creation of a spectra list of 1d spectra in fits format in a folder,
Correction/addition of the header.

Created on Fr Jul 19 14:48:03 2018
@author: lothar
"""

import glob
from astropy.io import fits


# Enter path and name of the spectra
files = input("Enter the path and name of the spectra: ")
filelist = glob.glob(files)

filelist.sort()

print("\nSpektrenliste: \n", filelist, "\n")
print("Anzahl der Spektren: ", len(filelist), "\n")

#   Header check first spectrum
hdul = fits.open(filelist[0])
print("\nInfo on the first spectrum of the list:")
hdul.info()
print("\nComplete header of the first spectrum of the list:")
print("\n", repr(fits.getheader(filelist[0], 0)))

Eingabe = input(
    "Would you like to change or add header data for all spectra in the series? \
If no enter 0, \
if header entry is a number: 1,  \
if header entry is a string: 2 : "
)

k = 0
KW = []
Wert = []
while Eingabe == '1' or Eingabe == '2':
    KW.append(input("Enter the header keyword: "))
    if Eingabe == "1":
        Wert.append(float(input("Enter the corresponding value: ")))
    elif Eingabe == "2":
        Wert.append(input("Enter the corresponding value: "))
    else:
        print("Incorrect input")
        break
    Eingabe = input(
        "Change other header entries? If no, enter 0, if yes, enter 1 or 2: "
    )
    k += 1


for i in range(len(filelist)):
    g = fits.open(filelist[i])
    for j in range(len(KW)):
        g[0].header[KW[j]] = Wert[j]
    fits.update(filelist[i], g[0].data, g[0].header)
    g.close()
