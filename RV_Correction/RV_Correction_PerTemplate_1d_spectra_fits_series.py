#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script corrects a series of 1d-spectra in fits format of an object
by an per Cross Correlation calculated RV [km/s].
Writes 1 fit-file of the RV-corrected spectra into the
working directory. The respective RV is noted in the the header of the
generated fit.

@author: Lothar Schanne
Stand 20211225
"""


from astropy.io import fits, ascii
import glob
from PyAstronomy import pyasl
import numpy as np
import matplotlib.pyplot as plt


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

# Alphabetical sorting. If the name is correct (e.g. 20180922-xyz.fit),
# this results in a temporal order.
filelist.sort()

# Print the list
print("\nList of spectra: \n")
print("Number of spectra: ", len(filelist), "\n")

# Select template
# Path and name of the template
tfile = input('Enter the path and file name of the template: ')

tf, theader = fits.getdata(tfile, header=True)
print(
    'Flux minimum and maximum in the template [ADU]: ', tf.min(), '  ', tf.max())
#   Check for the necessary header entries
print('Output of the header entries required for the wavelength calculation in the template:')
if 'NAXIS' in theader:
    print('Dimension, NAXIS:                        ', theader['NAXIS'])
else:
    print('This is not a 1d spectrum !')
if 'NAXIS1' in theader:
    tnax = theader['NAXIS1']
    print('Number of values (abscissa), NAXIS1:     ', tnax)
else:
    print('NAXIS1 missing in header !')
if 'CRVAL1' in theader:
    tcrval = theader['CRVAL1']
    print('Starting wavelength, CRVAL1:             ', tcrval)
else:
    print('CRVAL1 missing in header !')
if 'CDELT1' in theader:
    tcdel = theader['CDELT1']
    print('Wavelength step size, CDELT1:    ', tcdel)
else:
    print('CDELT1 missing in header !')
if 'CRPIX1' in theader:
    print('Reference pixel, CRPIX1: ', theader['CRPIX1'])
else:
    theader['CRPIX1'] = 1
#   Creating a numpy array with the wavelengths of the template
tw = np.ones(tnax, dtype=float)
tcrval = tcrval + (1 - theader['CRPIX1']) * tcdel
for i in range(tnax):
    tw[i] = tcrval + i * tcdel
print('Wavelength range of the template in angstroms: ', tw[0], tw[-1])
print('\n')

RV = np.zeros(len(filelist))


# Processing the filelist
for i in range(len(filelist)):
    #   Importing headers and data
    f, header = fits.getdata(filelist[i], header=True)
    print('\nFlux minimum and maximum in the target ' +
          filelist[i], f.min(), '  ', f.max())
    #   Check for the necessary header entries
    if 'NAXIS' in header:
        print('Dimension, NAXIS:                        ', header['NAXIS'])
    else:
        print('This is not a 1d spectrum !')
    if 'NAXIS1' in header:
        nax = header['NAXIS1']
        print('Number of values (abscissa), NAXIS1:     ', nax)
    else:
        print('NAXIS1 missing in header !')
    if 'CRVAL1' in header:
        crval = header['CRVAL1']
        print('Starting wavelength, CRVAL1:             ', crval)
    else:
        print('CRVAL1 missing in header  !')
    if 'CDELT1' in header:
        cdel = header['CDELT1']
        print('Wavelength step size, CDELT1:    ', cdel)
    else:
        print('CDELT1 missing in header  !')
    if 'CRPIX1' in header:
        print('Reference pixel, CRPIX1: ', header['CRPIX1'])
    else:
        header['CRPIX1'] = 1
    #   Creating a numpy array with the wavelengths of the target
    w = np.ones(nax, dtype=float)
    for j in range(nax):
        crval = crval + (1 - header['CRPIX1']) * cdel
        w[j] = crval + j * cdel
    print('Wavelength range in angstroms: ', w[0], w[-1])

    # Perform cross-correlation
    # The RV-range is (parameter 1) - to (parameter 2) km/s in
    # steps of (parameter 3) km/s.
    rv, cc = pyasl.crosscorrRV(
        w, f, tw, tf, -200., 200., 0.1, mode='doppler', skipedge=10)
    # Maximum of cross-correlation function
    maxind = np.argmax(cc)
    print("The cross-correlation function is maximum for dRV = ", rv[maxind],
          " km/s.")
    if rv[maxind] > 0.0:
        print("Red shift of the target compared to the template.")
    else:
        print("Blue shift of the target compared to the template")

    # Shift that spectrum
    flux_rv, wave_rv = pyasl.dopplerShift(
        w, f, -rv[maxind], edgeHandling="firstlast")

    header["CRVAL1"] = w[0]
    header['CRPIX1'] = 1

    newfile = filelist[i].rstrip(".fits") + "_RVcorrected_perTemplate.fit"
    header["RVCORR"] = (rv[maxind], "km/s, corrected")
    fits.writeto(newfile, flux_rv, header, overwrite=True,
                 output_verify="silentfix")

    RV[i] = rv[maxind]

    plt.figure(figsize=(10, 10))
    plt.title(filelist[i], size=10)
    plt.xlabel('Wavelength [Angström]')
    plt.ylabel('Flux normalized to the continuum')
    plt.plot(tw, tf, '-k', label='Template', linewidth=1.5)
    plt.plot(w, f, '-b', label='Original', linewidth=.5)
    plt.plot(w, flux_rv, '-r', label='adapted to template', linewidth=1.5)
    plt.legend(fontsize='x-small')
    plt.pause(1)
    # plt.savefig(filelist[i]+'RVcorr.pdf', format='pdf', format='pdf')
    plt.savefig(filelist[i]+'RVcorr.pdf', format='png', format='png')
    plt.close()


# Save the RVs as an ascii file
filename = filelist[i].rstrip(".fit")
ascii.write(
    [filelist, RV],
    filename + "_RVs_perTemplate.dat",
    overwrite=True,
    names=[
        "Spectrum",
        "RV"
    ],
    format="tab",
)
