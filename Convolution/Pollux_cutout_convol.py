#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Importing a synthetic spectrum in the form of a table (as available as .spec
in the Pollux database).
In the accompanying Pollux file Spektrum.txt the start wavelength and step
are specified.
Calculation of a selectable section (Panda's dataframe labeled 'range').

Plotted is the column NFLUX = normalized flux for the whole spectrum
and for range.
The range can be convolved with a standard deviation (stddev in the
unit of the step).
The convolved flux is stored together with the wavelength in a two-column
ascii table (spacer = tab, column names WAVE and FLUX).

Stand 20180823
@author: lothar
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve
from astropy.io import ascii

file = input('Enter path and file name: ')
lambda_min = float(input('Enter the start wavelength in the spectrum: '))
deltalambda = float(input('Enter the step width (pixels) in the spectrum: '))
std = float(input('Standard deviation for the convolution in multiples of\
 increment: '))

table = pd.read_fwf(file, names=['WAVE', 'AFLUX', 'NFLUX'], header=0)

wl_li = int(input('Wavelength range begin: '))
wl_re = int(input('Wavelength range end: '))
# bereich =  table.query('(wl_li<WAVE) & (WAVE<wl_re)')
ind_li = int((wl_li-lambda_min) / deltalambda)
ind_re = int((wl_re-lambda_min) / deltalambda)
bereich = table[ind_li:ind_re]

#   Convolve mit astropy.convolution
kernel = Gaussian1DKernel(stddev=std)
con = np.array(bereich['NFLUX'])
convoluted = convolve(con, kernel, normalize_kernel=True,
                      boundary='extend')

# Plots
fig, ax = plt.subplots(3)
fig.set_size_inches(15, 15)
plt.xlabel('Wavelength [Angström]')
fig.suptitle('Spectrum '+file)
plt.grid(True)
ax[0].plot(table['WAVE'], table['NFLUX'])
ax[1].plot(bereich['WAVE'], bereich['NFLUX'])
ax[2].plot(bereich['WAVE'], convoluted)
plt.savefig(file+'.png')
# plt.savefig(file+'.pdf')
plt.pause(.1)

# Saving the convolved
convol_file = pd.DataFrame(bereich['WAVE'], convoluted,
                           columns=['WAVE', 'NFLUX'])
name = file.rstrip('.spec')[0]+'_convolved.dat'
ascii.write([bereich['WAVE'], convoluted], name, overwrite=True,
            names=['WAVE', 'FLUX'], format='tab')
