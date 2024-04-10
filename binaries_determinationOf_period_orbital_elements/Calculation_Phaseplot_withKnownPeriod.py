#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculation_Phase_plot_withKnownPeriod.py

The script reads in an ascii file with the two columns (without column
heading !) JD (observation times) and RV (radial velocities).
The time points (JD) are convolved with the period to be entered (assumed to
be known)and the corresponding RV's are plotted against the phases.
The plot is saved as a png (or pdf).
The results are saved in a csv file with the column headings JD, RV and Phase.

20240218

@author: Lothar Schanne
"""

from astropy.io import ascii
import matplotlib.pyplot as plt
from PyAstronomy.pyasl import foldAt
import matplotlib.pylab as plt
import numpy as np
from astropy.io import ascii

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

filename = input("Enter the name of the data file: ")
ts = ascii.read(filename)
ts_JD = ts.columns[0]
ts_RV = ts.columns[1]

Periode = float(input('Enter the period in days: '))

phases = foldAt(ts_JD, Periode, ts_JD[0])

# Sort with respect to phase
# First, get the order of indices ...
sortIndi = np.argsort(phases)
# ... and, second, rearrange the arrays.
phases = phases[sortIndi]
RV = ts_RV[sortIndi]
JD = ts_JD[sortIndi]

fig = plt.figure()
plt.plot(phases*Periode, RV, 'ko', markersize=4)
plt.xlabel('Zeit [Tage]')
plt.ylabel('RV [km/s]')
plt.title('Phase plot ' + filename + ', calculated with period ' + str(Periode))
# plt.savefig(filename + '_Phaseplot_Periode_' + str(Periode), format='pdf')
plt.savefig(filename + '_Phaseplot_Periode_' + str(Periode), format='png')

plt.pause(.1)

# Saving the results as an ascii file (in csv format)
ascii.write(
    [JD, RV, phases],
    filename + "_Phase_Periode" + str(Periode) + "d.csv",
    overwrite=True,
    names=[
        "      JD     ",
        "      RV      ",
        '  Phase  ',
    ],
    format="csv",
)
