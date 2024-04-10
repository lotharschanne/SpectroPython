#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reading in a synthetic spectrum in the form of an ascii table with 2 columns
WAVE and FLUX overwritten.
Calculation of a selectable wavelength section. This range is then rebinned with
with a selectable step size and then additionally convoluted with a
selectable FWHM (apparatus profile).
The selected rebinned wavelength range is plotted and also the folded spectrum.
convolved spectrum.
The rebinned flux section of the original ascii file and the
rebinned and convolved flux is plotted. The data are saved in a
two-column ascii table (spacer = comma) and in one fits file.

Stand 20231113
@author: lothar
"""

import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import ascii, fits
from PyAstronomy.pyasl import binningx0dt


plt.switch_backend('Qt5Agg')  # Activate interactive graphics backend

file = input("PEnter path and file name: ")
table = ascii.read(file)

begin = float(
    input("Enter the start of the desired wavelength range in angstroms "))
end = float(
    input("Enter the end of the desired wavelength range: "))

fwhm = float(
    input("Enter the desired FWHM in angstroms: ")
)
step = float(
    input("Enter the desired increment in angstroms of the result spectrum: "))

newtable = table[table["WAVE"] >= begin]
newtable = newtable[newtable["WAVE"] <= end]

# Binning
data, d = binningx0dt(
    newtable["WAVE"],
    newtable["FLUX"],
    dt=step,
    x0=newtable["WAVE"].min(),
    removeEmpty=False,
)
Data = data.T
wave = Data[0]
flux = Data[1]

#   Convolve with astropy.convolution
kernel = Gaussian1DKernel(stddev=fwhm / 2.3 / step)
convoluted = convolve(flux, kernel, normalize_kernel=True, boundary="extend")

# Grafik
fig, ax = plt.subplots(2)
plt.xlabel("Wavelength [Angström]", fontsize=5)
fig.suptitle("Spectrum " + file, fontsize=5)
ax[0].plot(wave, flux, linewidth=0.2)
# ax[0].set_xlim(4000, 6700)
# ax[0].set_ylim(0.2, 1.1)
ax[0].tick_params(axis="both", labelsize=5)
ax[1].plot(wave, convoluted, linewidth=0.2)
ax[1].tick_params(axis="both", labelsize=5)
# ax[1].set_xlim(4000, 6700)
# ax[1].set_ylim(0.2, 1.1)
name = file.rstrip(".csv")  # anpassen !!!
plt.savefig(file + ".png")
plt.savefig(name + ".pdf")

# Save the trimmed ascii as dat and fits

# nicht convolviert
name = (
    file.rstrip(".csv") + "_" + str(int(begin)) +  # anpassen !!!!
    "_" + str(int(end)) + "_rebinned.dat"
)
ascii.write([wave, flux], name, overwrite=True,
            names=["WAVE", "FLUX"], format="tab")

name = (
    file.rstrip(".csv")  # anpassen !!!
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned.fits"
)

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CRVAL1"] = wave[0]
header["NAXIS1"] = len(wave)
header["CDELT1"] = step
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRPIX1"] = 1

fits.writeto(
    name, flux, header, overwrite=True, output_verify="silentfix",
)

# convolviert
name = (
    file.rstrip(".csv")   # abpassen !!!
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned_convolved.csv"
)

ascii.write(
    [wave, convoluted], name, overwrite=True, names=["WAVE", "FLUX"], format="csv"
)

name = (
    file.rstrip(".csv")   # anpassen
    + "_"
    + str(int(begin))
    + "_"
    + str(int(end))
    + "_rebinned_convolved.fits"
)

header = fits.Header()
header["SIMPLE"] = "T"
header["BITPIX"] = -32
header["NAXIS"] = 1
header["CRVAL1"] = wave[0]
header["NAXIS1"] = len(wave)
header["CDELT1"] = step
header["CUNIT1"] = "Angstrom"
header["CTYPE1"] = "Wavelength"
header["CRPIX1"] = 1

fits.writeto(
    name, convoluted, header, overwrite=True, output_verify="silentfix",
)
