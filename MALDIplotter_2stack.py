#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 13:23:27 2017

Takes two input mass spectra as text files and picks maxima.
Plots spectra vertically on different axes w/ shared x window & ticks.

@author: emeryusher
"""
# import python modules #
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

# import files into two separate arrays (index 0 contains m/z ratio; index 1
#   contains intensity. cannot contain headers)
array1 = np.loadtxt('spectrum1.txt')
array2 = np.loadtxt('spectrum2.txt')

# normalize data set to the maximum value in that dataset (this is especially good)
#   for comparing two separate samples

# np.max gives the maximum value in the array called "array1"
#   "array1[:,1]" means we are looking at the second column
#   in an array where the first (0) column is m/z and the second (1) column is intensity

array1[:,1] = array1[:,1]/ np.max(array1[:,1])
array2[:,1] = array2[:,1]/ np.max(array2[:,1])

# determine local extrema using the scipy feature we imported above
#   these values are stored as the array "max1" or "max2" for each dataset
# NB: there will be more than one maximum value given by this method; you can pick more maxima
#   by making the "order" number lower.
# google "agrelextrema" for more info on what's inside the parentheses
max1 = argrelextrema(array1[:,1], np.greater, order = 200)
max2 = argrelextrema(array2[:,1], np.greater, order = 200)


# plot intensity versus m/z using the "plt" module. The default is a line plot.
#   the syntax below results in an overlay of the two datasets.
# you can change the color of each line as you wish (google "matplotlib colors")
# change the label of each dataset; this feeds into the legend on the plot
fig, axs = plt.subplots(2, sharex = True)
fig.suptitle('mass spectra', fontsize = 14)
plt.setp(axs, xlim = (5000, 9000) , ylim = (-0.05, 1.5))
plt.setp(axs, ylabel = 'Intensity')
axs[0].plot(array2[:,0], array2[:,1], color = 'firebrick', linewidth = 1.0, label = 'sample 2')
axs[1].plot(array1[:,0], array1[:,1], color = 'teal', linewidth = 1.0, label = 'sample 1')

axs[0].legend(loc = 'best')
axs[1].legend(loc = 'best')

plt.xlabel('m/z')

# label the extrema for each data set within a specified window

# FOR all of the values in the "max1" array, IF a value is between [lower bound] and [upper bound],
#   then display the m/z value on the plot above the corresponding peak
for i in max1[0]:
    # change the less than/greater than values to fit the window you want to see (based on MW of your protein)
    if ((array1[:,0][i] > 7500 ) and (array1[:,0][i] < 8400 )) or ((array1[:,0][i] > 10000 ) and (array1[:,0][i] < 11000 )):
        axs[1].annotate(str(array1[:,0][i]), xy=(array1[:,0][i], array1[:,1][i]+0.02), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

for i in max2[0]:
    if ((array2[:,0][i] > 7500 ) and (array2[:,0][i] < 8400)):
        axs[0].annotate(str(array2[:,0][i]), xy=(array2[:,0][i], array2[:,1][i]+0.02), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

# uncomment the following line to save the stacked spectra as a png image
#plt.savefig('2mass_spectra.png', format='png', dpi=300)
