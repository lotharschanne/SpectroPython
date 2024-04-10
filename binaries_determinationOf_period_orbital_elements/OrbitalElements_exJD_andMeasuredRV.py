#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script reads an ascii file with two columns without headers containing the
times (JD, first column) and the RV's of a binary component (RV, second column).
The script calculates with the known (independently determined)
 parameters P and T the phase plot. For the remaining orbital elements, the
estimated values in lines 3 to 36 are generally sufficient.
The remaining orbital parameters K1, e, w1 and gamma are then optimized using
a curve fitting routine, so that the theoretical RVs can be calculated and
displayed in the phase plot for comparison with the measured RVs.
The optimized orbital parameters K1, e, w1, and gamma and their standard
deviations are printed out.

20231007
@author: lothar schanne
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Die RV-Tabelle im ascii-Format hat KEINE Spaltenüberschriften
table = input('Enter the file name with the RV data: ')
data = pd.read_table(table, names=['JD', 'RV'],
                     format='csv')  # Adapt format!

P = float(input('Enter the period to be used [in days]: '))

# # Initial orbital parameters of the binary system (estimated values)
# omega (w1) in degrees, gamma (system velocity in km/s) and K1 in km/s,
gamma = 0
K1 = 50
e = .5
w1 = 50
T = data.JD[0]


# Conversion of the time points relative to To and P (for phase plot):
JD = data.JD
zp = (JD - T) / P - (JD - T) // P

plt.plot(zp, data.RV, 'or', markersize=3., label='gemessene RV')


def pipe(zp, K1, e, w1, gamma):
    # Determination of the eccentric anomalies E by solution using Newton's method
    def newton(E):
        return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))

    E = np.ones(len(zp))
    v = np.ones(len(zp))
    vr1 = np.ones(len(zp))

    for n in range(len(zp)):
        rv = data.RV[n]
        while n >= 0:  # Iterative calculation of the eccentric anomaly E
            args1 = {'e': e, 'M': zp[n] * 2 * np.pi}
            args2 = {'e': e, 'M': zp[n] * 2 * np.pi}
            alt = E[n]
            E[n] = newton(alt)
            if abs(alt - E[n]) < 10e-12:
                break

        # Calculation of the true anomalies v:
        v[n] = 2 * np.arctan(np.tan(E[n] / 2) * np.sqrt(1 + e) / np.sqrt(1 - e))
        # Calculation of the radial velocity:
        vr1[n] = gamma + K1 * (e * np.cos(w1 * (np.pi/180))
                               + np.cos(w1 * (np.pi/180) + v[n]))
    return vr1


# Optimization of orbital parameters K1, e, w1 and gamma, please adjust bounds,
# if useful
popt, pcov = curve_fit(pipe, zp.to_numpy(),
                       data.RV.to_numpy(),
                       bounds=([0., 0., 0., -50.], [150., 1., 360., 50.]))

print('Optimized orbital parameters K1, e, w1 and gamma:\n', popt)
stdabw = np.sqrt(np.diag(pcov))
print('Standard deviations of the optimized parameters K1, e, w1 and gamma:\n',
      stdabw)


def pipe2(zp, popt):  # Calculation of the RVs with the optimized orbital parameters

    def newton(E):
        return E - (E - args1['e'] * np.sin(E) - args1['M']) / (1.0 - args2['e'] * np.cos(E))

    E = np.ones(len(zp))
    v = np.ones(len(zp))
    vr1 = np.ones(len(zp))

    for n in range(len(zp)):
        while n >= 0:  # Iterative calculation of the eccentric anomaly E
            args1 = {'e': popt[1], 'M': zp[n] * 2 * np.pi}
            args2 = {'e': popt[1], 'M': zp[n] * 2 * np.pi}
            alt = E[n]
            E[n] = newton(alt)
            if abs(alt - E[n]) < 10e-12:
                break

        # Calculation of the true anomalies v:
        v[n] = 2 * np.arctan(np.tan(E[n] / 2) *
                             np.sqrt(1 + popt[1]) / np.sqrt(1 - popt[1]))

        # Calculation of the radial velocity:
        vr1[n] = popt[3] + popt[0] * (popt[1] * np.cos(popt[2] * (np.pi/180))
                                      + np.cos(popt[2] * (np.pi/180) + v[n]))

    # Grafik:
    plt.plot(zp, vr1, 'b-',
             label='RV calculated with the optimized orbital parameters')
    plt.legend(fontsize='x-small')
    plt.xlabel('Phase')
    plt.ylabel('RV [km/s]')
    plt.title(table + ', P=' + str(P) + ' d' +
              '\nK=' + str(popt[0].round(1)) + '+-' +
              str(stdabw[0].round(1)) + ' km/s'
              ', e=' + str(popt[1].round(5)) + '+-' + str(stdabw[1].round(5))
              + ', omega=' + str(popt[2].round(1)) + '+-' +
              str(stdabw[2].round(1)) + '°'
              + ', gamma=' + str(popt[3].round(1)) + '+-' +
              str(stdabw[3].round(1)) + ' km/s',
              fontsize='small')
    plt.pause(1)
    grsave = input(
        'Would you like to save the graphic? Then answer with y :')
    if grsave == 'y':
        # plt.savefig(table+'_Phaseplot.pdf', format='pdf')
        plt.savefig(table+'_Phaseplot.png', format='png')
    plt.close('all')


z = np.linspace(0., 1., 1001)
pipe2(z, popt)
