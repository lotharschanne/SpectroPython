#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reads in a file (in csv or tab format) with the two columns JD
and the associated value of depending variable (e.g. RV's), but without column
headings.
Then performs a period analysis (PDM analysis) with both columns, displays the
periodogram and displays the period (in days).

Based on:
    https://pyastronomy.readthedocs.io/en/latest/pyTimingDoc/pyPDMDoc/pdm.html

Created on Thu Sep  7 16:00:55 2023

@author: lothar
"""


from astropy.io import ascii
import matplotlib.pylab as plt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import TimeSeries
from astropy import time
from astropy import units as u

plt.switch_backend('Qt5Agg')  # Interaktives GrafikBackend einschalten

filename = input("Enter the name of the data file: ")
ts = ascii.read(filename)
ts_JD = ts.columns[0]
ts_RV = ts.columns[1]

print('Enter the start and end of the period range, that is to be checked.')
pa = float(input('Enter the minimum P in the unit days: '))
pe = float(input('Enter the maximum P in the unit days: '))

# Get a ``scanner'', which defines the frequency interval to be checked.
# Alternatively, also periods could be used instead of frequency.
# Bitte anpassen !!!!!!!!!!!!!!!
S = pyPDM.Scanner(minVal=pa, maxVal=pe, dVal=0.001, mode="period")

# Carry out PDM analysis. Get frequency array
# (f, note that it is period, because the scanner's
# mode is ``period'') and associated Theta statistic (t).
P = pyPDM.PyPDM(ts_JD, ts_RV)
# Use 10 phase bins and 3 covers (= phase-shifted set of bins).
# Select the first parameter <= data points/2, also vary the second:
f1, t1 = P.pdmEquiBinCover(3, 3, S)

# Show the result
plt.figure(facecolor='white')
plt.title("Result of PDM analysis")
plt.xlabel("Frequency")
plt.ylabel("Theta")
plt.plot(f1, t1, 'bp-')
plt.pause(.1)


Periode = f1[t1.argmin()] * u.d
print('Periode = ', Periode.round(6))


# Create Time Series object
ts_JD = time.Time(ts.columns[0], format='jd')
serie = TimeSeries(time=ts_JD)
# Insert RV's
serie.add_column(ts_RV)

# Folding (Phasenplot)
a, b = P.phase(ts_JD, Periode)
fig = plt.figure()
plt.plot(a, ts_RV, 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nPDM phase diagram with the period ' +
          str(Periode.round(6))+' calculated')
plt.pause(.1)
# plt.savefig(filename + '_Phasediagram _PDM', format='pdf')
plt.savefig(filename + '_Phasediagram _PDM', format='png')
print('To exit the program, click on the last opened diagram (Figure 2)')
plt.waitforbuttonpress(-1)
plt.close('all')
