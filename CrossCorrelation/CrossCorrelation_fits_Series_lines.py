#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-correlation
derived from an example in PyAstronomy
https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/
crosscorr.html

A cross-correlation of a series of target spectra is performed with respect to a
template spectrum. Both are available as fits.
A spectrum section for the target spectra is queried (surrounding a line or a
center wavelength), which is used for the cross-correlation.
The RV and optionally the baryc. corrected RV are calculated. The data
are written to an ascii file.
The star and observer coordinates must be adjusted in the script.
Star coordinates can also be taken from the Internet.

Stand 20220408

@author: Lothar Schanne
"""


from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import glob
from astroplan import FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u

plt.switch_backend('Qt5Agg')  # Switch on interactive graphics backend


# ************* Coordinates of the observatory in the following form: *************
# longitude = 289.5967661 in Grad or [grad, min, sec]
# latitude = -24.62586583 in Grad or [grad, min, sec]
# altitude = 2635.43 in Meter

# ***************  PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# If this is omitted, the barycentric corrections are incorrect

# Coordinates of Berthold
# longitude = +7.4775
# latitude = 49.47527777778
altitude = 200
# Here alternative definition of the coordinates as a list in the form of degrees minutes seconds:
longitude = [7, 28, 39.0]
latitude = [49, 28, 31.0]


# Coordinates of  Wise Observatory in Israel
# longitude = +34.76333333
# latitude = 30.59583333
# altitude = 875

# Coordinates of Siegfried Hold
# longitude = +15.68461111
# latitude = 47.00161111
# altitude = 380


# Important !!!!!!!!!!!!!!!! :
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

# ***************  PLEASE CHANGE TO YOUR OWN COORDINATES ***************
# ************ Enter the coordinates of the star in the following form: *******
# If this is omitted, the barycentric corrections are incorrect
# ra2000 = 030.20313477 in Grad
# dec2000 = -12.87498346 in Grad
# or as a coordinate string RA DEC in the form "hh mmm ss +dd mm ss"

# # Coordinates of del Cep
# ra2000 = 337.29277083333335
# dec2000 = +58.415198
# coord = "22 29 10.265 +58 24 54.714"

# Coordinates of gam Cyg
# ra2000 = 305.55708
# dec2000 = +40.2566

# Coordinates of Betageuze
# ra2000 = 88.7929583
# dec2000 = 7.40705555

# Coordinates oftheta1 Ori C
# ra2000 = 83.81858333
# dec2000 = -5.389694444

# Coordinates of7 And
# ra2000 = 348.1375
# dec2000 = 49.40620275
# coord = "23 12 33 +49 24 22.3299"

# Coordinates of gam Cyg
# ra2000 = 305.55708
# dec2000 = 40.25666666

# Coordinates of OX Aurigae
# ra2000 = 103.25
# dec2000 = 38.86916
# coord = "06 53 01.41099 +38 52 08.9353"

# Coordinates of Polaris (alp UMi)
# coord = "02 31 49 +89 15 51"

# if type(coord) == str:
#     ra2000, dec2000 = pyasl.coordsSexaToDeg(coord)

# Reading in the star coordinates via the Internet
Frage = input(
    'Would you like to search for the star coordinates on the Internet? Then enter "y": ')
if Frage == 'y':
    star = input('Enter the name of the object star: ')
    ra2000 = FixedTarget.from_name(star).ra.value
    dec2000 = FixedTarget.from_name(star).dec.value
# ********************************************************************


# List of selectable lines. Please add.
Linien = {
    "O2 6883": 6883.818,
    "CIV 7726": 7726.2,
    "CIII 7037": 7037.25,
    "HeI 6678": 6678.149,
    "FeI 6678": 6677.9865,
    "FeI 6634": 6633.7492,
    "FeI 6609": 6609.1098,
    "H_alpha": 6562.817,
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

# Enter the line to be measured:
linie = input(
    "\nEnter the name of the line to be measured or a center wavelength \
        in angstroms: ")
if linie in Linien:
    wellenlaenge = Linien[linie]
else:
    wellenlaenge = float(linie)


# Enter the range around the line or center wavelength to be cut out
# and cross-correlated
bereich = float(
    input('Enter the range of the section around the line (in angstroms): '))

# Select template
# Path and name of the template
tfile = input("Enter the path and file name of the template: ")

#   Read header and data (flux from template in tf,
# header of the template stored in theader)
tf, theader = fits.getdata(tfile, header=True)
print(
    "\nMinimum and maximum in the template [ADU]: ", tf.min(), "  ", tf.max())

tnax = theader["NAXIS1"]
tcrval = theader["CRVAL1"]
step = theader["CDELT1"]
if 'CRPIX1' not in theader:
    theader['CRPIX1'] = 1

tw = np.ones(tnax, dtype=float)
tcrval = tcrval + (1 - theader["CRPIX1"]) * step
for i in range(tnax):
    tw[i] = tcrval + i * step

# Cut out range
index_wellenlaenge = int((wellenlaenge - tcrval) / step)
suchintervall = int((bereich + 10)/step)
t_intervallflux = np.zeros(suchintervall)
t_intervallwave = np.zeros(suchintervall)
for j in range(suchintervall):
    t_intervallflux[j] = tf[j + index_wellenlaenge - int(suchintervall / 2)]
    t_intervallwave[j] = (
        tcrval + (j + index_wellenlaenge - int(suchintervall / 2)) * step
    )

# Select the spectra (targets) to be correlated
# Path and name of the targets
# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the target spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results
# in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of target spectra:")
print(filelist)
print("\nNumber of target spectra: ", len(filelist), "\n")
frage_bary = input(
    'Would you like to correct the RV barycentrically? Then enter "y"')
print("Please wait. Calculation in progress")

RV = np.zeros(len(filelist))
RV_bc = np.zeros(len(filelist))
hjd = np.zeros(len(filelist))
corr = np.zeros(len(filelist))
obs_time = np.zeros(len(filelist))

for i in range(len(filelist)):
    f, header = fits.getdata(filelist[i], header=True)
    print(filelist[i])
    if frage_bary == 'j':
        if "JD-OBS" in header:
            obs_time[i] = float(header["JD-OBS"])
        elif "JD_OBS" in header:
            obs_time[i] = float(header["JD_OBS"])
        elif "MJD-OBS" in header:
            mjd = header["MJD-OBS"]
            obs_time[i] = mjd + 2400000.5
        elif "BAS_MJD" in header:
            obs_time[i] = float(header["BAS_MJD"])
        else:
            print("There is no observation time in the header of ",
                  filelist[i])
            break

        # Calculation of the heliocentric correction and time:
        corr[i], hjd[i] = pyasl.helcorr(
            longitude, latitude, altitude, ra2000, dec2000, obs_time[i],
            debug=False
        )

    # print('Minimum and Maximum im'+filelist[i]+' [ADU]: ', f.min(), '  ', f.max())
    nax = header["NAXIS1"]
    crval = header["CRVAL1"]
    step = header["CDELT1"]

    if 'CRPIX1' not in header:
        header['CRPIX1'] = 1

    #   Creating a numpy array with the wavelengths of the target
    w = np.ones(nax, dtype=float)
    crval = crval + (1 - header["CRPIX1"]) * step
    for j in range(nax):
        w[j] = crval + j * step

    # Cut out range
    lambda0 = float(header["CRVAL1"] + (1 - header["CRPIX1"]) * step)
    index_wellenlaenge = int((wellenlaenge - lambda0) / step)

    suchintervall = int(bereich/step)
    intervallflux = np.zeros(suchintervall)
    intervallwave = np.zeros(suchintervall)
    for j in range(suchintervall):
        intervallflux[j] = f[j + index_wellenlaenge - int(suchintervall / 2)]
        intervallwave[j] = (
            lambda0 + (j + index_wellenlaenge - int(suchintervall / 2)) * step
        )

    # Plot the interval
    fig = plt.figure(i)
    plt.title(filelist[i])
    plt.plot(intervallwave, intervallflux, '-r')
    plt.plot(t_intervallwave, t_intervallflux, '-b')
    plt.xlabel('Angström')
    plt.savefig(filelist[i]+'_cutout.png', format='png')
    plt.pause(1)

    # Perform cross-correlation.
    # The RV-range is (parameter 1) - to (parameter 2) km/s in
    # steps of (parameter 3) km/s.
    rv, cc = pyasl.crosscorrRV(
        intervallwave, intervallflux, t_intervallwave, t_intervallflux,
        -150, 150, 0.1, mode="doppler")

    # Maximum of cross-correlation function
    maxind = np.argmax(cc)
    RV[i] = rv[maxind]
    if frage_bary == 'y':
        RV_bc[i] = RV[i] + corr[i]
    else:
        RV_bc[i] = None

    # Plot CroCo-Funktion
    fig = plt.figure(2*i)
    plt.plot(rv, cc, "b-")
    plt.plot(rv[maxind], cc[maxind], "ro")
    plt.title("Cross-correlation function " + filelist[i])
    plt.grid(True)
    plt.xlabel("km/s")
    plt.savefig(filelist[i]+'_CroCo.png')
    plt.show()

print("Spectrum            ", '    JD     ',
      ' baryc. Corr. [km/s]', "     RV [km/s]", "        RV_bc [km/s]")
for i in range(len(filelist)):
    print(filelist[i], obs_time[i], corr[i], RV[i], RV_bc[i])

data = Table([filelist, obs_time, corr, RV, RV_bc], names=[
             "Spectrum", 'JD', 'baryz. Correction [km/s]', "RV [km/s]",
             'RV_bc [km/s]'])
ascii.write(data, "CC_RV_rangeWavelength_" +
            str(intervallwave[0].round()) + '_' +
            str(intervallwave[-1].round())
            + ".dat", overwrite=True, format="tab")

print('To exit the program, click on the last opened diagram')
plt.waitforbuttonpress(-1)
plt.close('all')
print("End of program")
