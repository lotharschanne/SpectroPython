# -*- coding: utf-8 -*-
"""
Reads in a file (in csv or tab format) with column headings. The
column headings for the independent variable (e.g. 'JD') and the dependent variable
variable (e.g. RV [km/s]) must be specified.
Then perform sa period analysis with both columns (Lomb-Scargle method and
string length method), graphically displays the pairs of values, the
periodogram and the phase plot and displays the period (in the unit of the
independent variable).


Lothar
20231229
"""


from PyAstronomy import pyTiming as pyt
import numpy as np
from astropy.timeseries import LombScargle, TimeSeries
from astropy.io import ascii
from astropy import units as u
from astropy import time
import matplotlib.pyplot as plt

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend

# Mean error of the RV measurements, please adjust
fehler = 3.  # km/s

filename = input("Enter the name of the data file: ")
x = input('Enter the name of the independent variable, e.g. "JD": ')
y = input('Enter the name of the dependent variable, e.g. "RV": ')

ts = ascii.read(filename)
ts_JD = time.Time(ts.columns[x], format='jd')
ts_RV = ts.columns[y]

print('Enter the start and end of the period range below, that is to be checked.')
pa = float(input('Enter the minimum P in the unit days: '))
pe = float(input('Enter the maximum P in the unit days: '))


# Time Series Object
serie = TimeSeries(time=ts_JD)

# RV's
serie.add_column(ts_RV)

fig = plt.figure()
plt.plot(serie.time.value, serie[y], "o")


# ************** Lomb-Scargle ex astropy.timeseries **************
frequency, power = LombScargle(serie.time.value, serie[y],
                               fehler).autopower(
    minimum_frequency=1/pe, maximum_frequency=1/pa)

fig = plt.figure()
plt.plot(frequency, power)
plt.title("Lomb-Scargle")
# plt.savefig('Lomb_Scargle', format='pdf')


power_max = power.max()
frequency_max = frequency[power.argmax()]
Periode = 1 / frequency_max * u.d
print("LombScargle method from astropy results in period = ", Periode.round(6))


# Folding (Phase)
serie_folded = serie.fold(period=Periode, epoch_time=ts_JD[0])
fig = plt.figure()
plt.plot(serie_folded.time.jd, serie_folded[y], 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nLomb-Scargle phase diagram with the period ' +
          str(Periode.round(5))+' calculated')
# plt.text(0, 0, 'Periode '+str(Periode.round(5)))
# plt.savefig(filename + '_Phaseplot_LombScargle', format='pdf')
plt.savefig(filename + '_Phaseplot_LombScargle', format='png')

plt.pause(.1)


# ****************** Further methods for period determination *************
jd = ts.columns[x]
rv = ts.columns[y]

# ************** String length method *****************
# Please adjust period range and number of steps
tps = (pa, pe, 10000)  # Tested period range [min, max, number of steps]
# Calculate string length
p, sl = pyt.stringlength_dat(jd, rv, tps)
per2 = p[np.argmin(sl)] * u.d
print('Period according to the string length method = ', per2.round(6))
# Show the string length.
fig = plt.figure()
plt.plot(p, sl, 'b.-')
plt.ylabel("String length")
plt.xlabel("Trial period")
plt.title('String length method')
plt.pause(.1)

# Folding (Phase)
serie_folded2 = serie.fold(period=per2, epoch_time=ts_JD[0])
fig = plt.figure()
plt.plot(serie_folded2.time.jd, serie_folded2[y], 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [km/s]')
plt.title(filename + '\nString length phase diagram with the period' +
          str(per2.round(5))+' calculated')
# plt.text(0, 0, 'Periode '+str(per2.round(5)))
plt.pause(.1)
# plt.savefig(filename + '_Phasenplot_Stringlenght', format='pdf')
plt.savefig(filename + '_Phasenplot_Stringlenght', format='png')

print('To exit the program, click on the last opened diagram (Figure 5)')
plt.waitforbuttonpress(-1)
plt.close('all')
