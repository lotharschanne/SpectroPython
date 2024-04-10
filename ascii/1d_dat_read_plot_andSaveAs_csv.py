#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d_dat_read_plot_csv_speichern.py

The script reads an ascii file. The ascii file has 2 columns
with the column headings 'WAVE' and 'FLUX'. The extension of the file is .dat.
The separator is determined automatically with guess=True.
Plots the file and saves as csv without column headings.

Created on Sun Jul 29 15:43:51 2018
@author: lothar
"""

import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import ascii


# Input of the asci file
# This one must have 2 separated columns WAVE and FLUX.
spectrum_name = input('Enter the path and name of the ascii file.: ')

# Read spectrum
spectrum = ascii.read(spectrum_name, guess=True)

# Create and plot spectrum as Pandas.series
pd_spectrum = pd.Series(spectrum['FLUX'], index=spectrum['WAVE'])
# print('click to determine the data points to normalize to the spectrum')
fig = plt.figure(figsize=(14, 10))
plt.plot(pd_spectrum)
plt.xlabel('Wavelength [Angström]')
plt.ylabel('ADU')
plt.title(spectrum_name)
plt.grid(True)

plt.show(block=True)

# Save as csv without column labeling
name = spectrum_name.strip('.dat')
name = name + '.csv'
pd_spectrum.to_csv(name)
