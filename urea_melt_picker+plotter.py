#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 14:46:39 2020

This script uses the formatted wavelength vs. intensity file for urea
denaturation data to pick the wavelength at which the emission intensity is
highest and plot the emission curve for each urea concentration on a single axis.

The output text file "max_wavs" feeds into

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from scipy.signal import argrelextrema

# import data from text file generated in "txt_file_organizers_v1.py"
xy_array = np.loadtxt('211123_WT_formatted.txt')

# text file containing ordered list of urea concentrations (in M)
urea_list = np.loadtxt('urea_list.txt')

urea_num = len(urea_list)

# split wavelength column from the rest of the array of spectra
wav = xy_array[:,0]
fluor = xy_array[:, 1:urea_num + 1]

# generate color map for temperature/urea concentration gradient
colormap = iter(cm.viridis(np.linspace(0, 1, urea_num)))

# plot all temperature spectra sequentially w/ a loop
for b in range (1, urea_num):
    plt.plot(wav, fluor[:,b], color = next(colormap), lw = 2)

# define color maps, normalize, and create key bar
cmap = cm.viridis
norm = mcolors.Normalize(vmin=min(urea_list), vmax=max(urea_list))
bar = cm.ScalarMappable(norm = norm, cmap = cmap)
bar.set_array(urea_list)

# axis stuff
plt.colorbar(bar, label = 'urea (M)')
plt.xlim(300, 400)
plt.xlabel('wavelength (nm)')
plt.ylabel('intensity (au)')
plt.title('MATH WT unfolding by urea' )

# uncomment the following line to save stacked urea melt plot
#plt.savefig('WT_stacked.png', format = 'png', dpi = 300)

first_max_index = argrelextrema(xy_array[:,1], np.greater, order = 500)

wavmax1 = np.transpose(fluor[first_max_index, :])
flat_max = wavmax1.ravel()

# this creates a list of the maximum intensity values
max_indices = []
max_Yarray = []
for i in range (1, urea_num+1):
    max_Yarray.append(xy_array[:,i].max())

# this finds the indices (max_indices) of the max intensity values
#   then turns that into an array of the values (temp)
#   finally, this values are stored as a list (max_Xarray)
for v in range(1, urea_num + 1):
    max_indices.append(argrelextrema(xy_array[:,v], np.greater, order = 500))
    temp = []
    for y in max_indices:
        temp.append(xy_array[:,0][y])

max_Xarray = []
for a in range(0, len(urea_list)):
    max_Xarray.append(temp[a][0])

# uncomment the following line to save the wavelengths of max. intensity
#np.savetxt('max_wavs.txt', max_Xarray, '%1.3f')

## plotting this stuff is optional (no fit, just the data points)
plt.figure()
plt.scatter(urea_list, max_Xarray, marker = 'o', color = 'crimson')
plt.xlabel('urea concentration (M)')
plt.ylabel('wavelength of maximum absorbance (nm)')
plt.title('WT unfolding by urea')

# uncomment the following line to save the plotted datapoints
#plt.savefig('WT_16pt_02nm.png', format = 'png', dpi = 600)
