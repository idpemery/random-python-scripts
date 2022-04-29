#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 13:23:27 2017

Takes five input mass spectra as text files and picks maxima.
Plots spectra vertically on different axes w/ shared x window & ticks.

@author: emeryusher
"""
#import python modules#

#numpy allows you to perform basic calculations, manipulate arrays, etc.
import numpy as np

#matplotlib is my favorite module for data presentation; there are others, but
#this one is the easiest, imo
import matplotlib.pyplot as plt

#from the scipy module, I am importing a feature that allows you to find extrema in a data set
from scipy.signal import argrelextrema

#import files into two separate arrays (one is unmodified, two is modified)

array1 = np.loadtxt('MR_ctrl_0_H1_1.txt')     #this spectrum will be on the bottom
array2 = np.loadtxt('MR_rxn1_0_H2_1.txt')
array3 = np.loadtxt('MR_rxn2_0_H3_1.txt')
array4 = np.loadtxt('MR_rxn3_0_H4_1.txt')
array5 = np.loadtxt('MR_rxn4_0_H5_1.txt')   #this spectrum will be on the top

# normalize data set to the maximum value in that dataset (this is especially good)
#   for comparing two separate samples

# np.max gives the maximum value in the array called "array1"
# "array1[:,1]" means we are looking at the second column
#   in an array where the first (0) column is m/z and the second (1) column is intensity

array1[:,1] = array1[:,1]/(np.amax(array1[:,1]))
array2[:,1] = array2[:,1]/(np.amax(array2[:,1]))
array3[:,1] = array3[:,1]/(np.amax(array3[:,1]))
array4[:,1] = array4[:,1]/(np.amax(array4[:,1]))
array5[:,1] = array5[:,1]/(np.amax(array5[:,1]))

# determine local extrema using the scipy feature we imported above

setord = 200

max1 = argrelextrema(array1[:,1], np.greater, order = setord)
max2 = argrelextrema(array2[:,1], np.greater, order = setord)
max3 = argrelextrema(array3[:,1], np.greater, order = setord)
max4 = argrelextrema(array4[:,1], np.greater, order = setord)
max5 = argrelextrema(array5[:,1], np.greater, order = setord)

# plot intensity versus m/z using the "plt" module. The default is a line plot.
# the syntax below results in an overlay of the two datasets.
# you can change the color of each line as you wish (google "matplotlib colors")
# change the label of each dataset; this feeds into the legend on the plot

fig, axs = plt.subplots(5, figsize = (5,10), sharex = True)
plt.setp(axs, xlim = (8800, 9800) , ylim = (-0.05, 1.5))
plt.setp(axs, ylabel = 'Intensity')

# don't judge me. sometimes working harder IS the smarter choice at the time
axs[0].plot(array5[:,0], array5[:,1], color = 'indigo', linewidth = 1.5, label = 'rxn 4')
axs[1].plot(array4[:,0], array4[:,1], color = 'forestgreen', linewidth = 1.5, label = 'rxn 3')
axs[2].plot(array3[:,0], array3[:,1], color = 'gold', linewidth = 1.5, label = 'rxn 2')
axs[3].plot(array2[:,0], array2[:,1], color = 'red', linewidth = 1.5, label = 'rxn 1')
axs[4].plot(array1[:,0], array1[:,1], color = 'k', linewidth = 1.5, label = 'apo')

axs[0].legend(loc = 'best')
axs[1].legend(loc = 'best')
axs[2].legend(loc = 'best')
axs[3].legend(loc = 'best')
axs[4].legend(loc = 'best')


plt.xlabel('m/z')

for ax in axs:
    ax.label_outer()

lbound = 9100
ubound = 9500

# label the extrema for each data set within a specified window

# FOR all of the values in the "max1" array, IF a value is between [lower bound] and [upper bound],
#   then display the m/z value on the plot above the corresponding peak
for i in max1[0]:
    # change the less than/greater than values to fit the window you want to see (based on MW of your protein)
    if ((array1[:,0][i] > lbound) and (array1[:,0][i] < ubound )):
        axs[4].annotate(str(array1[:,0][i]), xy=(array1[:,0][i], array1[:,1][i]+0.4), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

for i in max2[0]:
    if ((array2[:,0][i] > lbound) and (array2[:,0][i] < ubound)):
        axs[3].annotate(str(array2[:,0][i]), xy=(array2[:,0][i], array2[:,1][i]+0.4), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

for i in max3[0]:
    if ((array3[:,0][i] > lbound) and (array3[:,0][i] < ubound)):
        axs[2].annotate(str(array3[:,0][i]), xy=(array3[:,0][i], array3[:,1][i]+0.4), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

for i in max4[0]:
    if ((array4[:,0][i] > lbound) and (array4[:,0][i] < ubound)):
        axs[1].annotate(str(array4[:,0][i]), xy=(array4[:,0][i], array4[:,1][i]+0.4), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

for i in max5[0]:
    if ((array5[:,0][i] > lbound) and (array5[:,0][i] < ubound)):
        axs[0].annotate(str(array5[:,0][i]), xy=(array5[:,0][i], array5[:,1][i]+0.4), horizontalalignment = 'center', color = 'k', fontsize = 8, rotation=45)

# uncomment the following line to save the stacked spectra as a png image
#plt.savefig('5mass_spectra.png', format='png', dpi=300)
