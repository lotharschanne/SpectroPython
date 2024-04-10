#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script writes the observation date as JD at the beginning of the file name
of a series of spectra and saves the fits with the new name as fits.

20240201
@author: lothar schanne
"""
from astropy.io import ascii, fits
import glob


# Create file list. Spectra in a (sub)folder.
files = input("Path and name of the spectra (use wildcards) : ")
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
# print("\List of Spectra: \n")
# print(filelist, end='\n')
print("Number of spectra: ", len(filelist), "\n")


for i in range(len(filelist)):
    flux, header = fits.getdata(filelist[i], header=True)
    JD = header['JD']
    print('File name: ', filelist[i])
    name = 'JD_'+str(JD)+'_'+filelist[i]
    print('new file name: ', name)
    fits.writeto(name, flux, header, overwrite=True)
