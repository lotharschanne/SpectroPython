# -*- coding: utf-8 -*-
"""
Liest eine Datei (im csv- oder tab-Format) ein mit den beiden Spalten JD
und zugehöriger Wert (z.B. RV's), aber ohne Spaltenüberschriften.
Macht dann mit beiden Spalten eine Periodenanalyse (Lomb-Scargle-Methode),
zeigt grafisch die RV-Werte gegen die Zeit, das Periodogramm und den Phasenplot
und gibt die Periode (in der Einheit Tage) aus.

Lothar
20221107
"""


from PyAstronomy import pyTiming as pyt
import numpy as np
from scipy.signal import lombscargle as Lomb
from astropy.timeseries import LombScargle, TimeSeries
from astropy.io import ascii
from astropy import units as u
from astropy import time
import matplotlib.pyplot as plt

# Mittlerer Fehler der RV-Messungen, bitte anpassen
fehler = 3.  # km/s

filename = input("Geben Sie den Name des Datenfiles ein: ")
ts = ascii.read(filename)
ts_JD = time.Time(ts.columns[0], format='jd')
ts_RV = ts.columns[1]

# Time Series Objekt erstellen
serie = TimeSeries(time=ts_JD)

# RV's einfügen
serie.add_column(ts_RV)

fig = plt.figure()
plt.plot(serie.time.value, serie['col2'], "o")
# plt.savefig('Werteplot', format='pdf')

# Lomb-Scargle aus astropy.timeseries
# Alle Frequenzen (Einheit 1/Tage) im Skript bitte anpassen
frequency, power = LombScargle(serie.time.value, serie['col2'],
                               fehler).autopower(
    minimum_frequency=0.00001, maximum_frequency=.7
)

fig = plt.figure()
plt.plot(frequency, power)
plt.title("Lomb-Scargle")
# plt.savefig('Lomb_Scargle', format='pdf')


power_max = power.max()
frequency_max = frequency[power.argmax()]
Periode = 1 / frequency_max * u.d
print("LombScargle-Methode aus astropy ergibt Periode = ", Periode)


# Folding (Phase)
serie_folded = serie.fold(period=Periode, epoch_time=ts_JD[0])
fig = plt.figure()
plt.plot(serie_folded.time.jd, serie_folded['col2'], 'ko', markersize=3)
plt.xlabel('Time (days)')
plt.ylabel('RV [k/s]')
plt.title('Phasendiagramm mit der Periode '+str(Periode.round(5))+' berechnet')
# plt.text(0, 0, 'Periode '+str(Periode))
# plt.savefig('Phasenplot', format='pdf')

plt.show(block=True)


# weitere Methoden zur Periodenbestimmung
jd = ts.columns[0]
rv = ts.columns[1]

# Stringlängenmethode
# Periodenbereich und Schrittezahl bitte anpassen
tps = (13, 16, 100000)  # getesteter Periodenbereich [min, max, Schritteanzahl]
# Calculate string length
p, sl = pyt.stringlength_dat(jd, rv, tps)
per2 = p[np.argmin(sl)]
print('Periode nach der Stringlängenmethode = ', per2)
# Show the string length.
fig = plt.figure()
plt.plot(p, sl, 'b.-')
plt.ylabel("String length")
plt.xlabel("Trial period")
