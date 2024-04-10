#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BCalculates the heights of the objects in an object list (ascii file (csv) with
the object names object names, one in each line, see Objectlist.txt)
above the horizon for a planned point in time at an observation location and
plots them over a period of 24 hours. The graphics are saved as pdf in the
working directory.
The list of stars, the observer and the date/time must be entered.
The list of stars must exist in the working directory or in PYTHONPATH.

20221108

@author: lothar
"""


import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import EarthLocation
from astroplan import Observer
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astroplan.plots import plot_airmass, plot_altitude
from pytz import timezone
from astroplan import FixedTarget
import csv


# Data of the object
# coordinates = SkyCoord('9h23m55.42s', '+54d55m31.5s', frame='icrs')
# Read in the star coordinates via the Internet
objektliste = input('Enter the name of the list of object stars: ')

# gewünschtes Beobachtungsdatum
time = input('Enter the desired observation time.\
 In the form YYYY-MM-DD hh:mm:ss : ')
time = Time(time)

# Location of the observatory
ort = input('Enter the name of the observatory (observer): ')

# ++++++++++++++ List of observer locations ++++++++++++++
# please add analog if necessary
# Data of the upper observatory, here example in Israel, where NRES is installed
if ort == 'NRES':
    longitude = '34.763333'
    latitude = '30.595833'
    elevation = 1500 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='NRES Israel',
                   location=location,
                   timezone=timezone('Asia/Jerusalem'),
                   description="")

#  Innsbruck Observatorium in Anthing
if ort == 'Innsbruck':
    longitude = '13.008816'
    latitude = '47.886710'
    elevation = 800 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Flechas Innsbruck',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")

# Coordinates of Berthold (Glan-Münchweiler)
if ort == 'Berthold':
    longitude = '7.4775'
    latitude = '49.47527777778'
    elevation = 200 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Berthold Stober',
                        location=location,
                        timezone=timezone('Europe/Berlin'),
                        description="")

# Coordinates of Siegfried Hold
if ort == 'Siegfried':
    longitude = '5.68461111'
    latitude = '47.00161111'
    elevation = 380 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Siegfried Hold',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")

reader = csv.reader(open(objektliste))

for object in reader:
    object = object[0]
    coordinates = FixedTarget.from_name(object)
    # plot_airmass(coordinates, obs, time)
    fig = plt.figure()
    plot_altitude(coordinates, obs, time, airmass_yaxis=True)
    plt.title(object + ', ' + obs.name)
    fig.savefig(object+'.pdf')
    plt.pause(3)

print('To exit the program, click on the last opened diagram')
plt.waitforbuttonpress(-1)
plt.close('all')
