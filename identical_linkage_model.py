#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Identical Linkage Modeller
for competition FA experiment design

User must set the following:
    Kdy - dissociation constant of probe (fluorescently labeled binder) in M
    y - array of concentrations (in M) of probe
    Kdx - estimated dissociation constant for competitor molecule (in M)

Output is Kxobs in 1/M; maximize Kxobs while using high enough concentration
of probe molecule to achieve suitable S/N

@author: emeryusher
etusher@psu.edu
"""

# import required modules
import numpy as np
import matplotlib.pyplot as plt


## USER INPUT HERE ##
# Kdy is the dissociation constant of the macromolecule and fluorescent probe
#   molecule; must do direct titration to determine this
Kdy = 1e-6  # probe affinity for macromolecule in M

Kdx = 10e-6 # estimated competitor affinity for macromolecule in M

plot_title = 'plot_title'   # enter desired title for plot between apostrophes (no spaces)


## the rest of the script ##
Kay = 1 / Kdy

Kax = 1 / Kdx

# give a range of potential probe concentrations (in M)
# for reference, I get good S/N using 40 nM (40e-9 M) FITC-labeled probe
y = [10e-12, 10e-11, 10e-10, 10e-9, 10e-8, 10e-7, 10e-6, 10e-5, 10e-4, 10e-3]

# identical linkage model
Kxobs = Kax / (1 + Kay * np.array(y))

## plotting stuff ##
plt.scatter(y, Kxobs)
plt.xlim(10e-13, 0.1)
plt.xlabel('concentration of probe (M)')
plt.ylabel('Kax observed')
plt.xscale('log')
plt.title(str(plot_title))
plt.show()

saveif = input('Type "y" to save this figure as "' +str(plot_title)+ '.png". Type any other string to skip. ')

if saveif == "y":
    plt.savefig(str(plot_title)+ '.png',  format = 'png', dpi = 300)

else:
    print('Figure not save by user.')
