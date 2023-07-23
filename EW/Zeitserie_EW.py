#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet für eine Zeitserie im fits-Format die Äquivalentweite einer Linie.
Eingabe der Integrationsgrenzen grafisch-interaktiv oder manuell.
Die EW-Berechnung setzt voraus, dass für alle Spektren der Serie das gleiche
Wellenlängenintervall für die Berechnung des Integrals verwendet werden kann,
also keine wesentlichen RV-Änderungen stattfinden, oder wenn doch, dass die
Linie isoliert ist (also der Flux = 1 im Umfeld ist). Am besten baryzentrisch
korrigierte Spektren benutzen.

Created on Fri Oct 16 17:47:58 2020

@author: lothar
"""

import numpy as np
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import glob

plt.ion()

# Fileliste erstellen. Spektren in einem (Unter)Ordner.
files = input("Pfad und Name der Spektren (nutze wildcards) : ")
filelist = glob.glob(files)

# Alphabetisches Sortieren. Bei richtiger Namensgebung ergibt das eine
# zeitliche Ordnung
filelist.sort()

# Ausdruck der Liste
print("\nSpektrenliste: \n")
print(filelist)
print("Anzahl der Spektren: ", len(filelist), "\n")

EW_dict = {}

# Grafik erstes Spektrum
sp = fits.open(filelist[0], ignore_missing_end=True)
# print('\n\nHeader of the spectrum :\n\n', sp[0].header, '\n\n')

flux = np.array(sp[0].data)
wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

for i in np.arange(sp[0].header["NAXIS1"]):
    wave[i] = (
        sp[0].header["CRVAL1"]
        + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
    )
    # The list wave contains the wavelengths of the pixels.
# Close the fits-file:
sp.close()

# Plot the spectrum
# with plt.ion():
fig = plt.figure()
plt.plot(wave, flux)
plt.xlabel("Wavelength [Angström]")
plt.ylabel("ADU")
plt.title("Spectrum " + filelist[0])
plt.grid(True)
plt.show(block=False)

frage1 = input(
    "Möchten Sie die Integrationsgrenzen zahlenmäßig eingeben oder \
per Mausklick (grafisch)? Geben Sie m oder g ein: "
)

if frage1 == "g":
    # Interaktives Festlegen der Integrationsgrenzen
    # Falls die Darstellung des Spektrums interaktiv vergrößert wird die Eingabepunkte
    # jeweils mit der rechten Maustaste löschen, bis dann wirklich die linke
    # Seite der zu messenden Linie angeklickt wird.
    pts = []
    pts = np.asarray(plt.ginput(n=2, timeout=-1))
    plt.plot(pts[:, 0], pts[:, 1], "o", markersize=3)
    begin = pts[0, 0]
    end = pts[1, 0]
    print("Gewählter Integrationsbereich:", begin, " bis", end)
    print()

if frage1 == "m":
    # Eingabe der Integrationsgrenzen
    begin = float(
        input("Geben Sie die kurzwellige Wellenlänge-Integrationsgrenze ein: ")
    )
    end = float(
        input("Geben Sie die langwellige Wellenlänge-Integrationsgrenze ein: "))
    print()

EW = np.zeros(len(filelist))
JD = np.zeros(len(filelist))

# Abarbeiten der filelist
for k in np.arange(len(filelist)):
    sp = fits.open(filelist[k], ignore_missing_end=True)
    JD[k] = sp[0].header['JD']
    # Generation of arrays with the wavelengths and fluxes of the spectrum
    flux = np.array(sp[0].data)
    wave = np.ones(sp[0].header["NAXIS1"], dtype=float)

    for i in np.arange(sp[0].header["NAXIS1"]):
        wave[i] = (
            sp[0].header["CRVAL1"]
            + (i - sp[0].header["CRPIX1"] + 1) * sp[0].header["CDELT1"]
        )
    # The list wave contains the wavelengths of the pixels.
    # In the list flux the corresponding intensities.
    # Close the fits-file
    sp.close()

    for n in np.arange(len(wave)):
        if wave[n] <= begin and wave[n + 1] > begin:
            begin_n = n
            begin_wave = wave[n]
            begin_flux = np.median(flux[n - 2: n + 2])
            break
    for n in np.arange(len(wave)):
        if wave[n] <= end and wave[n + 1] > end:
            end_n = n
            end_wave = wave[n]
            end_flux = np.median(flux[n - 2: n + 2])
            break
    faktor = (end_flux - begin_flux) / (end_n - begin_n)
    ew = 0
    for m in np.arange((end_n - begin_n)):
        dif = begin_flux + m * faktor - flux[begin_n + m]
        ew += dif * sp[0].header["CDELT1"]
    EW[k] = ew

    print(filelist[k], "   EW =", EW[k])
    # Daten zum Weiterverarbeiten:
    EW_dict[filelist[k]] = (JD, begin_wave, end_wave,
                            begin_flux, end_flux, EW[k])

ascii.write(
    [filelist, JD, EW],
    'EW_' + str(int(begin)) + "_" + str(int(end)) + '.dat',
    names=["Spektrum", 'JD', "EW"],
    overwrite=True,
    format="tab",
)

fig = plt.figure()
plt.stem(JD, -EW)
plt.xlabel("JD")
plt.ylabel("-EW in Angstroem")
plt.title('EW ' + str(begin) + ' bis ' + str(end) + ' Angstroem')
plt.grid(True)
plt.savefig('EW_' + str(begin) + '_' + str(end) + '.pdf')
