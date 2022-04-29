#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulated FA titration (direct mode)

Uses Eq. 6 from Roehrl, Wang, & Wagner (2004) Biochemistry

User must set the following:
    Kd - dissociation constant of macromolecule and probe (in μM)
    M - macromolecule concnetrations for titration (>12 points; sample up to>10*Kd)
    probe - concentration of fluorescent probe molecule (uM)

@author: emeryusher
etusher@psu.edu
"""

import numpy as np
import matplotlib.pyplot as plt

# input dissociation constant (μM) and give array/file with macromolecule concentrations (μM)
Kd = 8
M = np.array([1000, 500, 250, 125, 62.5, 31.25, 15.625, 7.8125, 3.90625, 1.95312, 0.976562, 0.488281, 0.244141, 0.12207, 0.0610352, 0.0305176, 0.0152588, 0.00762939])
#M = np.loadtxt('sample_conc.txt')
# probe concentration (μM)
probe = 0.04

# model
FB =  0.14 * (0.03 + (((Kd + M + probe) - np.sqrt(((Kd + M + probe)**2) - 4 * M * probe)) / (2 * probe)))

# input some potential values for M to model the dynamic range
M1 = 1
M2 = 6
M3 = 12

# the 0.14 and 0.03 are scaling and y-translate factors, respectively, to convert FB to ~anisotropy
calcFB1 = 0.14 * (0.03 + (((Kd + M1 + probe) - np.sqrt(((Kd + M1 + probe)**2) - 4 * M1 * probe)) / (2 * probe)))
calcFB2 = 0.14 * (0.03 + (((Kd + M2 + probe) - np.sqrt(((Kd + M2 + probe)**2) - 4 * M2 * probe)) / (2 * probe)))
calcFB3 = 0.14 * (0.03 + (((Kd + M3 + probe) - np.sqrt(((Kd + M3 + probe)**2) - 4 * M3 * probe)) / (2 * probe)))


## PLOTTING ##
plt.figure(figsize = (4,4))

plt.plot(M, FB, color = 'grey', lw = 3)
plt.ylim(0, 0.16)
plt.xlim(0.0001, 2000)
plt.xscale('log')

plt.ylabel('anisotropy (au)')
plt.xlabel('macromolecule concentration (uM)')
plt.title('predicted FA range')

plt.scatter(M1, calcFB1, color = 'g', marker = '*', s = 250, label = str(M1)+' uM M')
plt.scatter(M2, calcFB2, color = 'm', marker = '*', s = 250, label = str(M2)+' uM M')
plt.scatter(M3, calcFB3, color = 'c', marker = '*', s = 250, label = str(M3)+' uM M')

plt.legend(loc = 'upper left')

a1 = round(calcFB1, 3)
a2 = round(calcFB2, 3)
a3 = round(calcFB3, 3)

plt.axhline(y = a1, linestyle = '--', color = 'g')
plt.axhline(y = a2, linestyle = '--', color = 'm')
plt.axhline(y = a3, linestyle = '--', color = 'c')

plt.tight_layout()
# uncomment the following line to save figure
#plt.savefig('predicted_FA_plot.png',  format = 'png', dpi = 300)
