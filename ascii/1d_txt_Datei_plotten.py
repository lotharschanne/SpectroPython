#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1d_txt_plotten.py

Das Skript liest eine ASCII-Datei ein. Die
Daten stehen in Spalten mit Überschriften.
Der Spaltenbezeichner der unabhängigen und abhängigen Variable werden abgefragt.
Die Daten lediglich geplottet.

Status 20180823
Author = Lothar Schanne
"""

import matplotlib.pyplot as plt
from astropy.io import ascii

name = input('Enter the path and name of the text file: ')
data = ascii.read(name, guess=True)

x = input('Geben Sie den Spaltennamen für die unabhängige Variable ein: ')
y = input('Geben Sie den Spaltennamen für die abhängige Variable ein: ')

# Plotting the spectrum
fig = plt.figure(figsize=(14, 10))
# plt.plot(data[x], data[y], 'bo')
plt.stem(data[x], -data[y])
plt.xlabel(x)
plt.ylabel(y)
plt.title(name)
plt.grid(True)

plt.savefig('EW, '+name+'.pdf')

# plt.show(block=True)
