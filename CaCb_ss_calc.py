#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 17:07:20 2021

This script calculates secondary structure using ΔδCα-ΔδCβ in reference
to random coil δCα and δCβ values.

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt

# text file must contain five columns in the following order (no headers):
#   residue # / experimental δCα / experimental δCβ / reference δCα / reference δCβ
data1 = np.loadtxt('cacb_data.txt')
residue = data1[:,0]

#calculate dCa and dCb
dCa = data1[:,1] - data1[:,3]
dCb = data1[:,2] - data2[:,4]

#calculate dCa - dCb
dCadCb = dCa - dCb

## setup the plot ##
ind = residue
width = 0.5
tick = np.arange(min(residue)+1, max(residue), 5)

fig = plt.figure(figsize = (6,2))
ax = fig.add_subplot(111)

apo = ax.bar(ind, dCadCb1, width, color = 'k')
phospho = ax.bar(ind+width, dCadCb1, width, color = 'magenta')

ax.set_ylabel('ΔδCα-ΔδCβ')
ax.set_title('secondary structure')
ax.set_xticks(tick)

# uncomment the following line to save the plot
#plt.savefig('ss_plot.png', format = 'png', dpi = 300)

plt.show()
