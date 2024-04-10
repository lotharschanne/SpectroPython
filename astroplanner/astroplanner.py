#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Berechnet die Höhe des Objektes über dem Horizont ab einem Zeitpunkt an einem
Beobachtungsort und plottet sie über einen Zeitraum von 24 Stunden.
Eingabe des Sterns, des Beobachters und des Datums nötig.

Created on Sat Feb  6 17:05:23 2021

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


frage = input('Would you like to enter the object celestial coordinates manually?\
 Then answer with y: ')

if frage == 'y':
    object = input('Enter the name of the object star: ')
    REC = input('Enter the rectascention in the form 9h23m55.42s: ')
    DEC = input('Enter the declination in the form +54d55m31.5s: ')
    coordinates = SkyCoord(REC, DEC, frame='icrs')
else:
    # Einlesen der Sternkoordinaten über das Internet
    object = input('Enter the name of the object star: ')
    coordinates = FixedTarget.from_name(object)

# Desired observation date
time = input('Enter the desired observation time.\
 In the form YYYY-MM-DD hh:mm:ss  : ')
time = Time(time)

# Location of the observatory
ort = input('Enter the name of the observatory (observer): ')

# ++++++++++++++ List of observer locations ++++++++++++++
# please add a new location analog if necessary
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

# Innsbruck Observatorium in Anthing
if ort == 'Innsbruck':
    longitude = '13.008816'
    latitude = '47.886710'
    elevation = 800 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Flechas Innsbruck',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")

# Koordinaten of Berthold (Glan-Münchweiler)
if ort == 'Berthold':
    longitude = '7.4775'
    latitude = '49.47527777778'
    elevation = 200 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Berthold Stober',
                        location=location,
                        timezone=timezone('Europe/Berlin'),
                        description="")

# Koordinaten of Siegfried Hold
if ort == 'Siegfried':
    longitude = '5.68461111'
    latitude = '47.00161111'
    elevation = 380 * u.m
    location = EarthLocation.from_geodetic(longitude, latitude, elevation)
    obs = Observer(name='Siegfried Hold',
                   location=location,
                   timezone=timezone('Europe/Berlin'),
                   description="")


# plot_airmass(coordinates, obs, time)
fig = plt.figure()
plot_altitude(coordinates, obs, time, airmass_yaxis=True)
plt.title(object + ', ' + obs.name)
grafik = object + '.pdf'
fig.savefig(grafik)
plt.pause(.1)

print('Zum beenden des Programms in das Diagramm klicken')
plt.waitforbuttonpress(-1)
plt.close('all')
