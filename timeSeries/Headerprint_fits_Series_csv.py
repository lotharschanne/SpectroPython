#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads the header data for a 1d spectra series in fits format and
writes them to separate ascii files (.csv).

Created on Fr Jul 19 14:48:03 2018
@author: Lothar Schanne
"""

import glob
from astropy.io import fits
import pandas as pd


files = input("Enter the path and name of the spectra: ")
filelist = glob.glob(files)

# Sortierung
filelist.sort()

# Ausdruck der Liste
print("\List of spectra: \n", filelist, "\n")
print("Number of spectra: ", len(filelist), "\n")

a = fits.open(filelist[0])
header = a[0].header
headerdict = dict(header)
file = filelist[0].rstrip(".fit")
pd_ser = pd.Series(headerdict.values(), index=headerdict.keys())
pd_fra = pd.DataFrame({file: pd_ser})

f = open(filelist[0].rstrip(".fit") + "_Header.csv", "w")
for i in headerdict:
    f.write(str(i) + "," + str(headerdict[i]) + "\n")
f.close()


for f in range(1, len(filelist)):
    g = fits.open(filelist[f])
    header = g[0].header
    headerdict = dict(header)
    f = open(filelist[f].rstrip(".fit") + "_Header.csv", "w")
    for i in headerdict:
        f.write(str(i) + "," + str(headerdict[i]) + "\n")
    f.close()
