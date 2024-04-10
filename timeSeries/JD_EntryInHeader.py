#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script reads in a series of fits-1d spectra and calculates the respective JD from
different header entries for the observation date and calculates the respective JD
and enters it in the header of the respective spectrum.

Created on Thu Dec  8 13:35:37 2022

@author: lothar
"""

import glob
from astropy.io import fits
from astropy.time import Time


# create filelist
# Enter path and name of the spectra
files = input("Enter the path and name of the spectra: ")
filelist = glob.glob(files)

filelist.sort()

# Printout of the list
# print("\nSpectrum list: \n", filelist, "\n")
print("Number of spectra: ", len(filelist), "\n")


# Processing the spectra list
for i in range(len(filelist)):
    g = fits.open(filelist[i])
    print()
    print(filelist[i], ' being processed')
    if 'JD' not in g[0].header:
        if 'OBSDATE' in g[0].header:
            try:
                t = Time(g[0].header['OBSDATE'])
                g[0].header['JD'] = t.jd
                print('OBSDATE found')
            except:
                print(filelist[i],
                      ': OBSDATE in wrong format')
        elif 'JD-OBS' in g[0].header:
            try:
                t = g[0].header['JD-OBS']
                g[0].header['JD'] = t
                print('JD-OBS found')
            except:
                print(filelist[i],
                      ': JD-OBS in wrong format')
        elif 'JD-MID' in g[0].header:
            try:
                t = g[0].header['JD-MID']
                g[0].header['JD'] = t
                print('JD-MID found')
            except:
                print(filelist[i],
                      ': JD-MID in wrong format')
        elif 'DATE-OBS' in g[0].header:
            try:
                ti = g[0].header['DATE-OBS']
                if ti.find('_') >= 0:
                    ti = ti.replace('_', ':')
                if ti.find(' ') >= 0:
                    ti = ti.replace(' ', '')
                t = Time(ti)
                g[0].header['JD'] = t.jd
                print('DATE-OBS found')
            except:
                print(filelist[i],
                      ': DATE_OBS in wrong format')

        elif 'DATE-BEG' in g[0].header:
            try:
                t = Time(g[0].header['DATE-BEG'])
                g[0].header['JD'] = t.jd
                print('DATE-BEG found')
            except:
                print(
                    filelist[i], ': DATE-BEG in wrong format')
        elif "MJD-OBS" in g[0].header:
            try:
                t = g[0].header["MJD-OBS"]
                g[0].header['JD'] = t + 2400000.5
                print('MJD-OBS found')
            except:
                print(
                    filelist[i], ': MJD-OBS in wrong format')
        elif "BAS-MJD" in g[0].header:
            try:
                t = g[0].header["BAS-MJD"]
                g[0].header['JD'] = t + 2400000.5
                print('BAS-MJD found')
            except:
                print(
                    filelist[i], ': BAS-MJD in wrong format')
        try:
            fits.update(filelist[i], g[0].data,
                        g[0].header, output_verify='ignore')
            print('JD entered in header')
        except:
            print('JD could not be entered in the header')
            pass
    else:
        print(filelist[i], ' Header entry JD already exists')
    g.close()
print('End of the program')
