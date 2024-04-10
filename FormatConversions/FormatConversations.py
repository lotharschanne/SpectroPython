#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Conversion JD -> Calendar date
Conversion calendar date -> JD
Conversion of hexagesimal RA and DEC to float
Conversion RA and DEC to float in hexagesimal RA and DEC

20221208

@author: lothar
"""

from PyAstronomy import pyasl
from astropy.time import Time

print("Selection of the task:")

typ = input(
    "If you want to convert a Julian date to calendar dates,\n\
 enter JD: \n\
If you want to convert a calendar date to a Julian date, enter\
 Enter CD: \n\
If you want to convert RA and DEC in hexagesimal form to float, enter \
RA: \n\
If you want to convert RA and DEC from float to hexagesimal form, enter \
DEC: \n\
"
)

if typ == "JD":
    # Convert JD to calendar date
    jd = float(input("Enter the JD: "))
    t = Time(jd, format='jd')
    print("Date: ", t.strftime('%H:%M:%S %d %b %Y'))
    print()

if typ == "CD":
    # Convert calendar date to JD
    kd = input(
        "Enter the calendar date in the form 2017-01-19T18:22:45 (=isot): \n")
    t = Time(kd, format='isot', scale='utc')
    print("Corresponding Julian date: ", t.jd)
    print("Corresponding reduced Julian date: ", t.mjd)
    print()

if typ == "RA":
    # Obtain decimal representation
    # The coordinate string. Valid formats are, e.g., “00 05 08.83239 +67 50 24.0135”
    # or “00:05:08.83239 -67:50:24.0135”.
    # Spaces or colons are allowed as separators for the individual components
    # of the coordinates.
    radec = input(
        "Enter RA and DEC in hexagesimal form. \n\
        Wie 00 05 08.83239 +67 50 24.01: \n"
    )
    print("Entered coordinates : ", radec)
    # Convert sexagesimal coordinates to decimal coordinates
    ra, dec = pyasl.coordsSexaToDeg(radec)
    print("Coordinates  [deg]: %010.6f  %+09.6f" % (ra, dec))

if typ == "DEC":
    ra = float(input("Enter RA as float: "))
    dec = float(input("Enter DEC as float: "))
    # Convert dezimal into sexagesimal representation
    print("RA, DEC", ra, dec)
    sexa = pyasl.coordsDegToSexa(ra, dec)
    print("Coordinates  [sexa]: ", sexa)
