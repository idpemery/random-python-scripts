#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:13:48 2021

This script uses the model described in:

 "Street TO, Courtemanche N, Barrick D 2008
 Protein Folding and Stability Using Denaturants.
 Biophysical Tools for Biologists,
 Volume One: In Vitro Techniques. pp. 295-325."

 (equation 9).

NB: this script optimizes deltaG and NOT m; the m value was fixed in this analysis
according to:

 "Myers JK, Pace CN, Scholtz JM.
 Denaturant m values and heat capacity changes: Relation to changes in
 accessible surface areas of protein unfolding.
 Protein Science. 1995; 4: 2138-2148."

using the crystallography change in accesible surface area upon unfolding
linear relation to m value.

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# provide [urea] and signal data files

urea = np.loadtxt('conc_list.txt')
raw_data = np.loadtxt('max_wavs.txt')

# compute "Yobs", a normalization from 0 (folded) to 1 (unfolded)
#   of the signal change that's being observed

Yobs = ((raw_data - np.min(raw_data)) / (np.max(raw_data) - np.min(raw_data)))

# constants (user sets temperature (in Kelvin) and m value (J/mol.M))
R = 8.314 #J/mol.K
T = 295 #K
m = 7700 #J/mol.M (where M is urea concentration in molar)

# define objective function to model the [urea]-dependent spectroscopic signal
    # G is free energy of unfolding in the absence of denaturant
    # Bd and Bn are m-values for denatured and native states (baseline slopes)
    # Yd and Yn are the y-intercepts for the denatured and native baselines

def denature(urea, G, Yd, Yn, Bd, Bn, C):

    exponent = (G + (m * urea)) / (R * T)
    num = Yd + (Bd * urea) + ((Yn + (Bn * urea)) * np.exp(-exponent))
    denom = 1 + np.exp(-exponent)

    Yobs = num / denom

    return C * Yobs

# user-inputted guesses for G, Yd, Yn, Bd, Bn, and C
    # suggested initial conditions:
    # G = -5000 J/mol
    # Yd = 1
    # Yn = 0
    # Bd = 0
    # Bn = 0
    # C = 1 (scaling factor; may omit)
guess = [-10000, 1, 0, 0, 0, 1]

# popt is where the fitted params are stored
popt, pcov = optimize.curve_fit(denature, urea, Yobs, guess)

# create 'simulated' x values to allow for curve smoothing
fitx = np.arange(np.min(urea), np.max(urea), 0.01)
# generate fit based on nonlinear least squares minimization (popt) params
fity = np.array(denature(fitx, *popt))

# determine error of fitted parameters (std. deviation)
perr = np.sqrt(np.diag(pcov))

# store G and m values for easy access (with conversion to kJ & rounding)
Ground = round(popt[0]) / 1000
Gstd = round(abs(perr[0])) / 1000

mkJ = m / 1000
color = 'dimgray'

## plotting stuff ##
plt.scatter(urea, Yobs, color = color)
plt.plot(fitx, fity, color = 'k', lw = 2)

plt.title('210928 WT urea denaturation')
plt.xlabel('urea concentration (M)')
plt.ylabel('relative signal (au)')
plt.xlim(-0.1, np.max(urea) + 0.1)
plt.ylim(-0.1, 1.1)

plt.annotate('$ΔG_{H2O}$ = ' +str(Ground)+ ' +/- ' +str(Gstd)+ ' kJ/mol', xy = (0.5, 0.5))
plt.annotate('m (fixed) = ' +str(mkJ)+ ' kJ/mol/M', xy = (0.5, 0.4))

## save the data ##
WT_data_out = data_out = np.column_stack([urea, Yobs])
WT_fit_out = np.column_stack([fitx, fity])

WT_params = [Ground, Gstd]

# uncomment the following lines to export the data, fit, and optimized parameters
#np.savetxt('rep2_data_out.txt', WT_data_out, '%1.3f')
#np.savetxt('rep2_fit_out.txt', WT_fit_out, '%1.3f')
#np.savetxt('rep2_params.txt', WT_params, '%1.3f')

# uncomment the following line to save the plot of the data w/ fit
#plt.savefig('WT_m7700.png',  format = 'png', dpi = 300)
