#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script is used to add additional Gaussian noise to a spectrum.
to a spectrum.

Created on Mon Dec  6 09:55:36 2021

@author: lothar
"""

from PyAstronomy import pyasl
import numpy as np
from astropy.io import fits, ascii
from random import gauss


file = input('Path and name of the spectrum (use wildcards): ')


flux, header = fits.getdata(file, header=True)
# The distribution of the noise is set with the sigma s of the Gaussian function
s = 0.2
for i in range(len(flux)):
    fehler = gauss(0, s)
    flux[i] = flux[i] + fehler


filename = file + '_noised_' + str(s)+'.fits'
fits.writeto(filename, flux,
             header, overwrite=True)
