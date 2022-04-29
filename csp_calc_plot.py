#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 16:07:00 2021

This script calculates Chemical Shift Perturbation (CSP) by determining the
difference between two points on a 2D plane. The following reference describes
the scaling constants "a_N" and "a_C" and considerations for their values:

"Williamson MP. Using chemical shift perturbation to characterise ligand binding.
Prog Nucl Magn Reson Spectrosc. 2013; 73: 1-16. PubMed PMID: 23962882"

This calculation can be easily adapted to 1H by removing the scaling factor.

@author: emeryusher
"""

#import modules
import numpy as np
import matplotlib.pyplot as plt

# load data (where state1 and state2 are the spectrum peaks to be compared)
#   column 1 must be the residue number, column 2 must be the N chemical shift,
#   and column 3 must be the C chemical shift

state1 = np.loadtxt('state1_shifts.txt')
state2 = np.loadtxt('state2_shifts.txt')

#specify scaling constants
a_N = 0.14
a_C = 0.30

#calculate delta N and delta C
N1 = state1[:,1]
N2 = state2[:,1]
deltaN = N2 - N1

C1 = state1[:,2]
C2 = state2[:,2]
deltaC = C2 - C1

# calculate CSP using distance function
CSP = np.sqrt((a_C * deltaC)**2 + (a_N * deltaN)**2)

## plotting stuff ##
plt.figure(figsize = (6,2))

plt.bar(state1[:,0], CSP, color = 'darkmagenta', edgecolor = 'white', linewidth = 1)
plt.plot(state1[:,0], CSP, color = 'darkmagenta', linewidth = 1)

plt.xlim(243, 283)
plt.ylabel('CSP (ppm)')
plt.xlabel('residue')
plt.title('CSPs upon phosphorylation')

# uncomment the following line to save the plot as an image
#plt.savefig('CSP_plot.png', format = 'png', dpi = 300)
