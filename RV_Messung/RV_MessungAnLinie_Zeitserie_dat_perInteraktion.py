#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest einen Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
tab-Format, 2 Spalten WAVE und FLUX). Festlegung des Linienminimums per
Interaktion und Bestimmung der Radialgeschwindigkeit aus dem Minimum.
Plottet alle Spektren zum markieren von bis zu 2 Linien-Minima (im Falle eines
Doppelsterns SB2) und gibt die ermittelten Daten (RV und Apex) als ascii-Dateien
(Komma-separiert, als .csv) aus. Die RV's sind nicht baryzentrisch korrigiert.'

Stand 20221105

@author: lothar schanne
"""

import numpy as np
from astropy.io import ascii
import glob
from PyAstronomy import pyasl
import matplotlib.pyplot as plt

# plt.style.use('seaborn-whitegrid')


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
    "CIV_7726": 7726.2,
    "CIII_7037": 7037.25,
    "HeI_6678": 6678.149,
    "FeI_6678": 6677.9865,
    "FeI_6634": 6633.7492,
    "FeI_6609": 6609.1098,
    "H_alpha": 6562.817,
    "FeI_6546": 6546.24,
    "FeII_6516": 6516.0783,
    "FeI_6463": 6462.725,
    "FeII_6456": 6456.3805,
    "FeI_6417": 6416.9386,
    "FeI_6412": 6411.6592,
    "FeI_6408": 6408.0272,
    "FeI_6400": 6400.0008,
    "FeI_6394": 6393.6009,
    "SiII_6371": 6371.36,
    "SiII_6347": 6347.10,
    "FeI_6265": 6265.1335,
    "FeI_6256": 6256.3611,
    "FeI_6255": 6254.581,
    "FeI_6253": 6252.555,
    "FeII_6248": 6247.559,
    "FeI_6233": 6232.6408,
    "FeI_6231": 6230.7226,
    "FeI_6213": 6213.299,
    "FeI_6200": 6200.3125,
    "FeI_6192": 6191.558,
    "FeI_6180": 6180.2038,
    "NiI_6177": 6176.81,
    "NiI_6175": 6175.367,
    "FeI_6173": 6173.3352,
    "FeII_6170": 6169.816,
    "FeI_6170": 6169.597,
    "FeI_6164": 6163.5441,
    "CaI_6162": 6162.17,
    "FeII_6149": 6149.231,
    "FeII_6148": 6147.734,
    "HgII_6142": 6141.773,
    "FeI_6142": 6141.7316,
    "FeI_6137": 6137.286,
    "CaI_6122": 6122.22,
    "NiI_6108": 6108.12,
    "CaI_6103": 6102.72,
    "FeI_6065": 6065.482,
    "FeI_6056": 6056.0043,
    "FeI_6027": 6027.0505,
    "FeI_6024": 6024.0576,
    "FeI_6020": 6020.1688,
    "CuII_6013": 6013.411,
    "FeIII_5920": 5920.0,
    "CII_5920": 5919.6,
    "FeIII_5919": 5918.960,
    "Na_D1": 5895.924,
    "Na_D2": 5889.951,
    "HeI_5876": 5875.989,
    "FeI_5860": 5859.608,
    "FeI_5862": 5862.357,
    "HI_5861": 5861.35,
    "CIV_5812": 5812.140,
    "CIV_5802": 5801.510,
    "CIII_5696": 5696.000,
    "OIII_5592": 5592.370,
    "OIII_5593": 5592.37,
    "FeI_5576": 5576.0884,
    "FeI_5573": 5572.842,
    "FeI_5570": 5569.6177,
    "FeI_5567": 5567.3907,
    "FeI_5566": 5565.7036,
    "FeI_5560": 5560.2112,
    "FeI_5558": 5557.9818,
    "FeI_5555": 5554.8947,
    "FeI_5526": 5525.5439,
    "FeI_5498": 5497.5157,
    "TiII_5491": 5490.7,
    "FeI_5447": 5446.8743,
    "FeI_5430": 5429.6964,
    "TiII_5419": 5418.8,
    "FeI_5415": 5415.1989,
    "FeII_5412": 5411.970,
    "HeII_5411": 5411.524,
    "FeI_5406": 5405.7749,
    "TiI_5404": 5404.11,
    "FeI_5383": 5383.3688,
    "TiII_5381": 5381.03,
    "FeI_5367": 5367.466,
    "FeII_5363": 5362.9698,
    "FeI_5307": 5307.36,
    "FeI_5302": 5302.299,
    "TiII_5262": 5262.14,
    "FeI_5233": 5232.94,
    "MgI_5183": 5183.6042,
    "MgI_5172": 5172.6843,
    "MgI_5167": 5167.3216,
    "FeI_5162": 5162.2725,
    "FeII_5154": 5154.242,
    "FeI_5097": 5096.9977,
    "FeI_5075": 5074.748,
    "TiII_5072": 5072.25,
    "FeI_5065": 5065.0181,
    "VI_5062": 5061.79,
    "FeI_5018": 5018.4354,
    "HeI_5016": 5015.6783,
    "FeI_5007": 5007.275,
    "FeI_5002": 5001.8633,
    "HeI_4922": 4921.929,
    "H_beta": 4861.332,
    "HeI_4713": 4713.143,
    "HgII_4687": 4686.563,
    "HeI_4686": 4685.682,
    "CIII_4650": 4650.16,
    "CIII_4647": 4647.400,
    "MgI_4571": 4571.0956,
    "FeII_4542": 4541.985,
    "HeI_4542": 4541.59,
    "He_4452": 4541.59,
    "HeI_4471": 4471.688,
    "HeI_4388": 4387.928,
    "CIII_4388": 4388.24,
    "H_gamma": 4340.468,
    "NIII_4200": 4200.020,
    "HI_4102": 4101.737,
    "SiIV_4089": 4088.863,
    "HeI_4026": 4026.191,
}
print("Linienauswahl: ", Linien.keys())


# Eingabe der zu messenden Linie:
linie = input("Geben Sie den Namen der zu messenden Linie ein: ")
wellenlaenge = Linien[linie]

# Abarbeiten der filelist, Einlesen von flux und Wellenlängen, Auswahl des Flux um die Linie:

# Definition von Variablen
RV1 = np.zeros(len(filelist))
RV2 = np.zeros(len(filelist))
apex1 = np.zeros(len(filelist))
apex2 = np.zeros(len(filelist))
linien_minwave1 = np.zeros(len(filelist))
miniwave1 = np.zeros(len(filelist))
miniwave2 = np.zeros(len(filelist))
miniflux1 = np.zeros(len(filelist))
miniflux2 = np.zeros(len(filelist))


for i in range(len(filelist)):
    ts = ascii.read(filelist[i], format="tab")

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall = 4
    # ************************************************************

    intervall = np.extract(
        abs(ts.columns[0] - wellenlaenge) < suchintervall, ts)
    a = np.zeros(len(intervall))
    b = np.zeros(len(intervall))
    for j in range(len(intervall)):
        a[j], b[j] = intervall[j]

    fig = plt.figure(i)
    plt.plot(a, b, "-", linewidth=0.5)
    plt.title(filelist[i].rstrip(".dat") + "_" + linie)

    print("\nSpektrum ", filelist[i])
    print("Klicke das Minimum der ersten Linie an: ")
    pts = np.asarray(plt.ginput(n=1, timeout=-1))
    RV1[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299792
    apex1[i] = pts[0, 1]
    print("Erste Linie, RV =", RV1[i], "Apex =", apex1[i])

    frage = input(
        "Möchten Sie das Minimum einer zweiten Linie anklicken? Dann y eingeben: "
    )
    if frage == "y":
        print("Klicke die zweite Linie an")
        pts = np.asarray(plt.ginput(n=1, timeout=-1))
        RV2[i] = (pts[0, 0] - wellenlaenge) / wellenlaenge * 299792
        apex2[i] = pts[0, 1]
        print("Zweite Linie, RV =", RV2[i], "Apex =", apex2[i])
    else:
        RV2[i] = np.NaN
        apex2[i] = np.NaN

    # plt.savefig(filelist[i].rstrip(".dat") + "_" + linie + ".png")


# Abspeichern als ascii-Datei
ascii.write(
    [filelist, RV1, apex1, RV2, apex2],
    linie + "_RV_interaktiv" + ".csv",
    overwrite=True,
    names=["Spektrum", "RV1", "Apex1", "RV2", "Apex2"],
    format="csv",
)
