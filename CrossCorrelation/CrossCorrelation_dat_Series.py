#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cross-correlation
derived from an example in PyAstronomy
https://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/
crosscorr.html

A cross-correlation of a series of target spectra is performed with respect to a
template spectrum is performed. Both are available as ascii files, format .dat.
The calculated RVs are printed out and saved in a file.

Stand 20220408

@author: Lothar Schanne
"""


from PyAstronomy import pyasl
import numpy as np
import matplotlib.pylab as plt
from astropy.io import ascii
from astropy.table import Table
import glob
from linetools.spectra.xspectrum1d import XSpectrum1D


obj = input("Enter the heading for the graphics (without blank):")

# Select template
# Path and name of the template
template_name = input("Enter the path and file name of the template: ")
template = np.loadtxt(template_name, skiprows=1)
print("Minimum and maximum in the template [ADU]: ",
      template.min(), "  ", template.max())

w = template[:, 0]
f = template[:, 1]

dt = (w[-1] - w[0]) / len(w)

template_newdata, dt_template = pyasl.binningx0dt(
    w, f, dt=dt, x0=w[0], useBinCenter=True
)

# print("Template begin und end: ", template_newdata[0, 0], template_newdata[-1, 0])
print("Template begin and end: ", w[0], w[-1])

# Select spectra to be correlated (target)
# Path and name of the targets
# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the target spectra (use wildcards) : ")
filelist = glob.glob(files)

# Sort alphabetically. If the spectrum files are named correctly, this results # in a temporal order.
filelist.sort()

# Printout of the list for control purposes.
print("\nList of target spectra:")
print(filelist)
print("\nNumber of target spectra: ", len(filelist), "\n")
print("Please wait. Calculation in progress")

RV = np.zeros(len(filelist))
jd = np.zeros(len(filelist))

for i in range(len(filelist)):
    target = np.loadtxt(filelist[i], skiprows=1)
    # print('Minimum and Maximum im'+filelist[i]+' [ADU]: ', f.min(), '  ', f.max())

    #   Generate one numpy array each with the wavelengths and flux of the target
    tw = target[:, 0]
    tf = target[:, 1]
    binnumber = int((tw[-1] - tw[0]) / dt_template)

    target_newdata, dt_target = pyasl.binningx0dt(
        tw, tf, x0=tw[0], nbins=binnumber, useBinCenter=True
    )

    print(filelist[i])
    print("Target Beginn und Ende:",
          target_newdata[0, 0], target_newdata[-1, 0])

    # Perform cross-correlation.
    # The RV-range is (parameter 1) - to (parameter 2) km/s in
    # steps of (parameter 3) km/s.
    rv, cc = pyasl.crosscorrRV(
        target_newdata[:, 0],
        target_newdata[:, 1],
        w,
        f,
        -200.0,
        200.0,
        0.1,
        mode="doppler",
        skipedge=100,
    )

    # Maximum of cross-correlation function
    maxind = np.argmax(cc)
    RV[i] = rv[maxind]
#
    # Plot CroCo-Funktion
    fig = plt.figure()
    plt.plot(rv, cc, "b-")
    plt.plot(rv[maxind], cc[maxind], "ro")
    plt.title(obj + " : Cross-correlation function " + filelist[i])
    plt.grid(True)
    plt.xlabel("km/s")
    plt.savefig(filelist[i]+'_CroCo.png')
    plt.show()

print("Spectrum", 'JD', "RV [km/s]")
for i in range(len(filelist)):
    print(filelist[i], jd[i], RV[i])

data = Table([filelist, jd, RV], names=["Spectrum", 'JD', "RV"])
ascii.write(data, "RV_List_" + obj + ".dat", overwrite=True, format="tab")

print('To exit the program, click on the last opened diagram')
plt.waitforbuttonpress(-1)
plt.close('all')

print("End of Program")
