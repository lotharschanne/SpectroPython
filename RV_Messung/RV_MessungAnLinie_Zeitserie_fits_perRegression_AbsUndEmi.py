#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Liest einen Spektrenkatalog ein (Zeitserie von normierten 1d-Spektren im
fits-Format). Berechnet aus dem Beobachtungszeitpunkt und den (anzupassenden)
Koordinaten des Beobachters und Objekts die heliozentrische Korrektur,
fitted die gewählte Linie im Minimum(Absorption)- oder Maximum(Emission)bereich 
per Regression und bestimmt aus dem heliozentrisch korrigierten Minimum/Maximum 
die heliozentrisch korrigierte Radialgeschwindigkeit RV.
Plottet und speichert alle fittings und gibt die ermittelten Daten als ascii-Datei
(Komma-separiert, als .csv) aus.

Stand 20221105

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import glob
import matplotlib.pylab as plt
from PyAstronomy import pyasl
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

# plt.style.use('seaborn-whitegrid')


# ************* Koordinaten des Observatory in folgender Form: *************
# longitude = 289.5967661 in Grad oder [grad, min, sec]
# latitude = -24.62586583 in Grad oder [grad, min, sec]
# altitude = 2635.43 in Meter

# ***************  BITTE AUF EIGENE KOORDINATEN ÄNDERN ***************
# Falls das versäumt wird sind die baryzentrischen Korrekturen falsch

# Koordinaten vom Berthold
# longitude = +7.4775
# latitude = 49.47527777778
# altitude = 200
# longitude = [7, 28, 39.0]
# latitude = [49, 28, 31.0]


# Koordinaten des Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Koordinaten von Siegfried Hold
longitude = +15.68461111
latitude = 47.00161111
altitude = 380


# Wichtig !!!!!!!!!!!!!!!! :
# esign : int, optional, {-1,0,1}
# Explicit sign with -1 representing negative sign, +1 representing positive
# sign, and 0 indicating no explicit sign specification. The explicit sign is
# necessary if negative southern coordinates are specified but d is 0 and,
# thus, cannot carry the sign.

if type(longitude) == list:
    longitude = pyasl.dmsToDeg(longitude[0], longitude[1], longitude[2])

if type(latitude) == list:
    latitude = pyasl.dmsToDeg(latitude[0], latitude[1], latitude[2], esign=+1)

# ********************************************************************

# ***************  BITTE AUF EIGENE KOORDINATEN ÄNDERN ***************
# ************ Eingabe der Koordinaten des Sterns in folgender Form: *******
# Falls das versäumt wird sind die baryzentrischen Korrekturen falsch
# ra2000 = 030.20313477 in Grad
# dec2000 = -12.87498346 in Grad
# oder als Koordinaten-String RA DEC in der Form "hh mmm ss +dd mm ss"

# # Koordinaten von del Cep
# ra2000 = 337.29277083333335
# dec2000 = +58.415198
# coord = "22 29 10.265 +58 24 54.714"

# Koordinaten von gam Cyg
# ra2000 = 305.55708
# dec2000 = +40.2566

# Koordinaten von Betageuze
# ra2000 = 88.7929583
# dec2000 = 7.40705555

# Koordinaten von theta1 Ori C
# ra2000 = 83.81858333
# dec2000 = -5.389694444

# Koordinaten von 7 And
# ra2000 = 348.1375
# dec2000 = 49.40620275
# coord = "23 12 33 +49 24 22.3299"

# Koordinaten von gam Cyg
# ra2000 = 305.55708
# dec2000 = 40.25666666

# Koordinaten von OX Aurigae
# ra2000 = 103.25
# dec2000 = 38.86916
# coord = "06 53 01.41099 +38 52 08.9353"

# Koordinaten von Polaris (alp UMi)
# coord = "02 31 49 +89 15 51"

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Einlesen der Sternkorrdinaten über das Internet
Frage = input(
    'Möchten Sie die Sternkoordinaten im Internet suchen lassen? Dann "y" eingeben: ')
if Frage == 'y':
    star = input('Geben Sie den Namen des Objektsterns ein: ')
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************

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
    "FeIII 5920": 5920.0,
    "CII 5920": 5919.6,
    "FeIII 5919": 5918.960,
    "Na D1": 5895.924,
    "Na D2": 5889.951,
    "HeI_5876": 5875.989,
    "FeI 5860": 5859.608,
    "FeI 5862": 5862.357,
    "HI 5861": 5861.35,
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
    "Geben Sie den Namen der zu messenden Linie oder eine Wellenlänge in\
    Angström ein: ")
if linie in Linien:
    wellenlaenge = Linien[linie]
else:
    wellenlaenge = float(linie)

frage_minmax = input('Geben Sie e ein, wenn es sich um eine Emissionslinie \
handelt, andernfalls a: ')

# Eingabe der Systemgeschwindigkeit
systemgeschwindigkeit = float(
    input("Geben Sie eine Systemgeschwindigkeit in km/s  ein: ")
)
systemgeschwindigkeit = systemgeschwindigkeit / 299772 * wellenlaenge

frage_bary = input(
    'Möchten Sie die RV baryzentrisch korrigieren? Dann Eingabe von "j"')

# Abarbeiten der filelist, Einlesen von flux und header, Auswahl des Flux um die Linie:

# Definition von Variablen
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))
RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
miniflux = np.zeros(len(filelist))

#################################################################
# Regression, Grad anpassen (2 oder 4 oder 6)
grad = 4
#################################################################

for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    if frage_bary == 'j':
        if "JD-OBS" in header:
            obs_time[i] = float(header["JD-OBS"])
        elif 'JD' in header:
            obs_time[i] = float(header["JD"])
        elif "JD_OBS" in header:
            obs_time[i] = float(header["JD_OBS"])
        elif "MJD-OBS" in header:
            mjd = header["MJD-OBS"]
            obs_time[i] = mjd + 2400000.5
        elif "BAS_MJD" in header:
            obs_time[i] = float(header["BAS_MJD"])
        else:
            print("Es ist kein Beobachtungszeitpunkt im Header von ",
                  filelist[i])
            break
        # Berechnung der heliozentrischen Korrektur und Zeit:
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i],
            debug=False
        )
        print("\n" + filelist[i] + ":")
        print("Beobachtungszeitpunkt: ", obs_time[i])
        print("Barycentric correction [km/s]: ", corr[i].round(2))
    else:
        corr[i] = 0

    if "CRPIX1" in header:
        refpix = int(header["CRPIX1"])
    else:
        refpix = 1

    step = float(header["CDELT1"])

    lambda0 = float(header["CRVAL1"]) - (refpix - 1) * step

    index_wellenlaenge = int(
        (wellenlaenge - lambda0 + systemgeschwindigkeit) / step)

    # **********************************************************
    # Breite des Suchintervalls, Breite in Angström, anpassen.
    suchintervall = int(10 / step)
    # ************************************************************
    intervallflux = np.zeros(suchintervall)
    intervallwave = np.zeros(suchintervall)
    for j in range(suchintervall):
        intervallflux[j] = flux[j +
                                index_wellenlaenge - int(suchintervall / 2)]
        intervallwave[j] = (
            lambda0 + (j + index_wellenlaenge - int(suchintervall / 2)) * step
        )

    if frage_minmax == 'a':
        linienminimum_wave_index = intervallflux.argmin()
    else:
        linienminimum_wave_index = intervallflux.argmax()

    ##########################################################################
    # Das folgende Intervall anpassen, bei verrauschten Linien breit (-20,+21)
    # bei nicht verrauschten, schmalen Linien klein (-2,+3)
    ##########################################################################
    linienminimum_intervall = np.arange(
        linienminimum_wave_index - 5, linienminimum_wave_index + 6
    )

    # Wellenlängen und Flux direkt um das Linienminimum
    wl = intervallwave[linienminimum_intervall]
    fl = intervallflux[linienminimum_intervall]

    model = np.poly1d(np.polyfit(wl, fl, grad))

    polyline = np.linspace(wl[0], wl[-1], 100)
    modflux = model(polyline)

    # Linienextremum berechnen
    if frage_minmax == 'a':
        miniflux[i] = modflux.min()
        miniwave = polyline[modflux.argmin()]
    else:
        miniflux[i] = modflux.max()
        miniwave = polyline[modflux.argmax()]
    miniwave_bc = miniwave * (1 + corr[i] / 299792)

    # Plotten
    fig = plt.figure()
    plt.plot(intervallwave, intervallflux, "o-")
    plt.plot(polyline, modflux)
    plt.plot(miniwave, miniflux[i], "o", color="black")
    plt.title(filelist[i])
    plt.xlabel("Wellenlänge [Angström]")
    plt.ylabel("relative Intensität")
    plt.savefig(filelist[i].rstrip(".fit") + "_" + linie + "_Regression.png")

    # RV ohne baryzentrische Korrektur
    RV[i] = (miniwave - wellenlaenge) / wellenlaenge * 299792
    print("nicht korrigierte RV: ", RV[i].round(2))
    # RV bc-korrigiert, sysv nicht berücksichtigt:
    if frage_bary == 'j':
        RV_bc[i] = (miniwave_bc - wellenlaenge) / wellenlaenge * 299792
        print("baryzentrisch korrigierte RV: ", RV_bc[i].round(2))
    else:
        RV_bc[i] = None


# # Plot der RV's
# fig = plt.figure()
# plt.plot(obs_time, RV, 'bo', markersize=1)

# Abspeichern als ascii-Datei
ascii.write(
    [filelist, obs_time, corr, RV, RV_bc, miniflux],
    linie + "_RV_Regression_grad" + str(grad) + ".csv",
    overwrite=True,
    names=[
        "Spektrum",
        "JD",
        "BC [km/s]",
        "RV [km/s]",
        "RV_bc [km/s]",
        "Flux im Extremum",
    ],
    format="csv",
)
