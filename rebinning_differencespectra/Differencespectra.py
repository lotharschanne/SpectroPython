#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Differenzspektren.py

Reads tab spectra of a series with the columns "WAVE" and "FLUX and forms
differential spectra to the specified reference spectrum, which are then saved
as a tab file.
All spectra used must have the same wavelength range and the same step size.
All difference spectra are displayed graphically in one plot each, which can
also be saved.

20231115

@author: lothar
"""
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import glob


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the dat spectra (use wildcards): ")
filelist = glob.glob(files)

# Alphabetical sorting. With the correct naming, this results in a
# chronological order
filelist.sort()

# Printout of the list
print("\nSpectrum list: \n")
print(filelist, end='\n')
print("Number of spectra: ", len(filelist), "\n")

mean_name = input("Name of the dat file for the reference spectrum: ")
mean_spec = ascii.read(mean_name, format="tab")
print("Reference spectrum used", mean_name)

for i in range(len(filelist)):
    spec = ascii.read(filelist[i], format="tab")
    diff = np.zeros(len(spec["WAVE"]))
    for j in range(len(spec["WAVE"])):
        diff[j] = spec["FLUX"][j] - mean_spec["FLUX"][j]
    filename = filelist[i].rsplit(".")[0] + "_diffspectrum" + ".dat"
    ascii.write(
        [spec["WAVE"], diff[:]],
        filename,
        overwrite=True,
        names=["WAVE", "FLUX"],
        format="tab",
    )
    fig = plt.figure(figsize=(10, 10))
    plt.title(filename, fontsize=10)
    plt.xlabel("Wavelength [Angstroem]")
    plt.ylabel("Differenz Flux")
    plt.ylim(-.3, .45)  # anpassen !!!!!!!!!!!!!!!!!!
    plt.grid(True)
    plt.plot(spec["WAVE"], diff)
    plt.pause(.2)
    fig.savefig("Differenzspectrum" +
                filelist[i].rsplit('.')[0] + ".png", format='png')
    plt.close()
