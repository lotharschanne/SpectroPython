#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script calculates from the orbital parameters assumed to be known and
a series of time points (in JD, are read in from an ascii table)
the theoretical radial velocities of the 2 components.
Writes the results to a csv-file named RVs.csv and plots the phase diagram.

20231007
@author: lothar schanne
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

datenherkunft = input('Would you like to import the times from an\
ascii file? If yes, then enter y: ')
if datenherkunft == 'y':
    t = input('Enter the file name of the time points: ')
    t = np.loadtxt(t)
else:
    # Example times of the observations (JD), !!!!! adjust !!!!!
    t = np.array([
        2.459971531251342967e+06,
        2.459972730320548173e+06,
        2.459977545927382074e+06,
        2.459990486079244874e+06,
        2.459994471415292937e+06,
        2.459995491419071797e+06])

# # Orbital parameters of the binary system,
# omegas (w) in Grad),
# P und T in JD, gamma (System speed) and K1 und K2 in km/s,
# please adapt to the specific case #############################
OP = {
    'gamma': 0.82407167,
    'K1': 31.8565906,
    'K2': 29.8379765,
    'T': 2459685.699575991,
    'e': 0.00212,
    'w1': 0.5168984*180/np.pi,
    'w2': (0.5168984 + np.pi)*180/np.pi,
    'P': 71.656017}

# Conversion of the times relative to To and P (in phase):
zp = (t - OP['T']) / OP['P'] -\
    (t - OP['T']) // OP['P']

# Mean anomalies at the points in time zp:
M = zp * 2 * np.pi

# Determination of the eccentric anomalies E by solution using Newton's method
def newton(E):
    return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))


E = np.ones(len(M))
v = np.ones(len(M))
vr1 = np.ones(len(M))
vr2 = np.ones(len(M))


for n in range(len(M)):
    iter = 0
    while n >= 0:  # Iterative calculation of the eccentric anomaly E
        args1 = {'e': OP['e'], 'M': M[n]}
        args2 = {'e': OP['e'], 'M': M[n]}
        alt = E[n]
        E[n] = newton(alt)
        iter += 1
        # print(n, 'Iter.: ', iter, 'E :', E[n])
        if abs(alt - E[n]) < 10e-10:
            break

    # Calculation of the true anomalies v:
    v[n] = 2 * np.arctan(np.tan(E[n] / 2) * np.sqrt(1 + OP['e'])
                         / np.sqrt(1 - OP['e']))
    # Calculation of the two radial velocities:
    vr1[n] = OP['gamma'] + OP['K1'] *\
        (OP['e']*np.cos(OP['w1']*(np.pi/180))
         + np.cos(OP['w1']*(np.pi/180) + v[n]))

    vr2[n] = OP['gamma'] + OP['K2']\
        * (OP['e']*np.cos(OP['w2']*(np.pi/180))
           + np.cos(OP['w2']*(np.pi/180) + v[n]))
#    Result printout
for i in range(len(vr1)):
    print('JD:', t[i], 'rv1:', vr1[i].round(2), 'rv2:', vr2[i].round(2))

ascii.write([t, vr1, vr2], 'RVs.csv', overwrite=True, names=['JD', 'RV1', 'RV2'],
            format="csv")

# Grafik:
plt.plot(zp, vr1, 'or', label='vr1')
plt.plot(zp, vr2, 'ob', label='vr2')
plt.legend()
plt.pause(1)
grsave = input('Would you like to save the graphic ? Then answer with y :')
if grsave == 'y':
    # plt.savefig('Phaseplot.pdf')
    plt.savefig('Phaseplot.png')
plt.close('all')
