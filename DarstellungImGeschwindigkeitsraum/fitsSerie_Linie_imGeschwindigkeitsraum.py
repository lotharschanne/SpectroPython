#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Einlesen einer Serie von heliozentrisch korrigierten Spektren im fits-Format,
Darstellung einer wählbaren Linie im Geschwindigkeitsraum in zwei Plots:
    Erstens alle Spektren übereinander geplottet.
    Zweitens alle Spektren mit einem wählbaren offset übereinander geplottet.
Die plots können als pdf abgespeichert werden.

20230216

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt


# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print("Anzahl der Spektren: ", len(filelist), "\n")


# Liste der wählbaren Linien:
Linien = {
    "O2 6883": 6883.818,
    "CIV 7726": 7726.2,
    "CIII 7037": 7037.25,
    'FeI 6837': 6837.0057,
    "HeI 6678": 6678.149,
    "FeI 6678": 6677.9865,
    "FeI 6634": 6633.7492,
    "FeI 6609": 6609.1098,
    "H alpha": 6562.817,
    "FeI 6546": 6546.24,
    "FeII 6516": 6516.0783,
    "FeI 6463": 6462.725,
    "FeII 6456": 6456.3805,
    "FeI 6417": 6416.9386,
    "FeI 6412": 6411.6592,
    "FeI 6408": 6408.0272,
    "FeI 6400": 6400.0008,
    "FeI 6394": 6393.6009,
    "SiII 6371": 6371.36,
    "SiII 6347": 6347.10,
    "FeI 6265": 6265.1335,
    "FeI 6256": 6256.3611,
    "FeI 6255": 6254.581,
    "FeI 6253": 6252.555,
    "FeII 6248": 6247.559,
    "FeI 6233": 6232.6408,
    "FeI 6231": 6230.7226,
    "FeI 6213": 6213.299,
    "FeI 6200": 6200.3125,
    "FeI 6192": 6191.558,
    "FeI 6180": 6180.2038,
    "NiI 6177": 6176.81,
    "NiI 6175": 6175.367,
    "FeI 6173": 6173.3352,
    "FeII 6170": 6169.816,
    "FeI 6170": 6169.597,
    "FeI 6164": 6163.5441,
    "CaI 6162": 6162.17,
    "FeII 6149": 6149.231,
    "FeII 6148": 6147.734,
    "HgII 6142": 6141.773,
    "FeI 6142": 6141.7316,
    "FeI 6137": 6137.286,
    "CaI 6122": 6122.22,
    "NiI 6108": 6108.12,
    "CaI 6103": 6102.72,
    "FeI 6065": 6065.482,
    "FeI 6056": 6056.0043,
    "FeI 6027": 6027.0505,
    "FeI 6024": 6024.0576,
    "FeI 6020": 6020.1688,
    "CuII 6013": 6013.411,
    'SiII 5979': 5978.93,
    "FeIII 5920": 5920.0,
    "CII 5920": 5919.6,
    "FeIII 5919": 5918.960,
    "Na D1": 5895.924,
    "Na D2": 5889.951,
    "HeI 5876": 5875.989,
    "FeI 5860": 5859.608,
    "FeI 5862": 5862.357,
    "CIV 5812": 5812.140,
    "CIV 5802": 5801.510,
    "CIII 5696": 5696.000,
    "OIII 5592": 5592.370,
    "OIII 5593": 5592.37,
    "FeI 5576": 5576.0884,
    "FeI 5573": 5572.842,
    "FeI 5570": 5569.6177,
    "FeI 5567": 5567.3907,
    "FeI 5566": 5565.7036,
    "FeI 5560": 5560.2112,
    "FeI 5558": 5557.9818,
    "FeI 5555": 5554.8947,
    "FeI 5526": 5525.5439,
    "FeI 5498": 5497.5157,
    "TiII 5491": 5490.7,
    "FeI 5447": 5446.8743,
    "FeI 5430": 5429.6964,
    "TiII 5419": 5418.8,
    "FeI 5415": 5415.1989,
    "FeII 5412": 5411.970,
    "HeII 5411": 5411.524,
    "FeI 5406": 5405.7749,
    "TiI 5404": 5404.11,
    "FeI 5383": 5383.3688,
    "TiII 5381": 5381.03,
    "FeI 5367": 5367.466,
    "FeII 5363": 5362.9698,
    "FeI 5307": 5307.36,
    "FeI 5302": 5302.299,
    "TiII 5262": 5262.14,
    "FeI 5233": 5232.94,
    "MgI 5183": 5183.6042,
    "MgI 5172": 5172.6843,
    "MgI 5167": 5167.3216,
    "FeI 5162": 5162.2725,
    "FeII 5154": 5154.242,
    "FeI 5097": 5096.9977,
    "FeI 5075": 5074.748,
    "TiII 5072": 5072.25,
    "FeI 5065": 5065.0181,
    "VI 5062": 5061.79,
    "FeI 5018": 5018.4354,
    "HeI 5016": 5015.6783,
    "FeI 5007": 5007.275,
    "FeI 5002": 5001.8633,
    "HeI 4922": 4921.929,
    "H beta": 4861.332,
    "HeI 4713": 4713.143,
    "HgII 4687": 4686.563,
    "HeI 4686": 4685.682,
    "CIII 4650": 4650.16,
    "CIII 4647": 4647.400,
    "MgI 4571": 4571.0956,
    "FeII 4542": 4541.985,
    "HeI 4542": 4541.59,
    "He 4452": 4541.59,
    "HeI 4471": 4471.688,
    "HeI 4388": 4387.928,
    "CIII 4388": 4388.24,
    "H gamma": 4340.468,
    "NIII 4200": 4200.020,
    "HI 4102": 4101.737,
    "SiIV 4089": 4088.863,
    "HeI 4026": 4026.191,
}
print("Linienauswahl: ", Linien.keys())

# Eingabe der zu messenden Linie:
linie = input(
    "Geben Sie den Namen der zu messenden Linie oder eine Wellenlänge in \
Angström ein: ")
if linie in Linien:
    wellenlaenge = Linien[linie]
else:
    wellenlaenge = float(linie)

bereich = float(input('Geben Sie den darzustellenden Geschwindigkeitsbereich um \
die gewählte Linie in km/s ein: '))

# Eingabe der Systemgeschwindigkeit:
systemgeschwindigkeit = float(
    input("Geben Sie die Systemgeschwindigkeit in km/s ein: "))
systemgeschwindigkeit = systemgeschwindigkeit / 299772 * wellenlaenge

# Parameter für den Abstand zwischen den Spektren im offset Plot fig1:
offset = float(input("Please enter the desired offset: "))
obj = input("Please enter the object name: ")


# Abarbeiten der Spektrenliste:
for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1
    step = float(header["CDELT1"])
    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step
    JD = header['JD']

    wave = np.zeros(header['NAXIS1'])
    for j in range(header['NAXIS1']):
        wave[j] = lambda0 + j * step

    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

    vstep = step / wellenlaenge * 299772
    pixbereich = int(bereich / vstep)

    fluxbereich = flux[(index_wellenlaenge - pixbereich)
                        :(index_wellenlaenge + pixbereich)]
    wellenlaengenbereich = wave[(
        index_wellenlaenge - pixbereich):(index_wellenlaenge + pixbereich)]
    geschwindigkeitsbereich = (wellenlaengenbereich - wellenlaenge) / \
        wellenlaenge * 299772

    plt.figure(1, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich, linewidth=.5)
    plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
    plt.ylabel('relative Intensität')
    plt.title(obj + ', ' + linie)

    plt.figure(2, figsize=(7, 10))
    plt.plot(geschwindigkeitsbereich, fluxbereich +
             i * offset, "-", label=JD, linewidth=1)
    # Beschriftung der einzelnen Spektren:
    plt.text(geschwindigkeitsbereich[1], 1.0 +
             i * offset, format(JD, '.2f'), ha="left", size=7)
    plt.xlabel('Geschwindigkeit relativ zur Laborwellenlänge in km/s')
    plt.ylabel('relative Intensität')
    plt.title(obj + ', ' + linie)
    plt.axvline(x=0, color='k', linewidth=.05)
    plt.grid(visible=True, axis='x')
    plt.pause(.01)

frage = input(
    'Möchten Sie die Grafiken abspeichern als pdf? Wenn ja "y" eingeben: ')
if frage == 'y':
    plt.figure(1)
    plt.savefig('Overplot_' + obj + '_'+linie + '.pdf')
    plt.savefig('Overplot_' + obj + '_'+linie + '.png')
    plt.figure(2)
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '.pdf')
    plt.savefig('OverplotMitOffset' + obj + '_'+linie + '.png')

plt.figure(1).clear()
plt.figure(2).clear()
