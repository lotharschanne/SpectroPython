#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script reads a csv file with the columns JD and RV (without column headings).
Using the period and T0 parameters specified in the script (please adjust)
the observation times are then convolved, the phase plot is displayed and
the results are saved in an ascii file (format tab, columns 'JD', 'Phase', 'RV').

Created on Tue Nov 14 19:10:56 2023

@author: lothar
"""
from PyAstronomy import pyasl
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

Periode = 2.867328  # Please adjust
T0 = 2441773.49     # Please adjust

filename = input("Enter the name of the data file: ")
ts = ascii.read(filename)

jd = ts.columns[0]
rv = ts.columns[1]


# Folding (Phase)
phases = pyasl.foldAt(jd, Periode, T0)
fig = plt.figure()
plt.plot(phases, rv, 'ko', markersize=3)
plt.xlabel('Phase')
plt.ylabel('RV [km/s]')
plt.title(filename.rsplit('.csv')[0] + 'Phase diagram with the period ' +
          str(Periode)+' d caculated')
# plt.text(0, 0, 'Periode '+str(Periode.round(5)))
# plt.savefig(filename + '_Phaseplot', format='pdf')
plt.savefig(filename + '_Phaseplot', format='png')

# plt.show(block=True)
plt.pause(.1)

data = Table([jd, phases, rv], names=['JD', 'Phase', 'RV'])
ascii.write(data, filename.rsplit('.csv')[0] + "_Phasis.dat",
            overwrite=True, format="tab")
