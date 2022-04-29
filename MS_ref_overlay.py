#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 21 13:23:27 2017

Plots two mass spectra on one set of axes and notates m/z of the peak(s)
with the highest relative intensity.

@author: emeryusher
"""

# import python modules for plotting etc.
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema

# import files into two separate arrays (one is the control; one is the reaction)
temp1 = np.loadtxt('chloro_ref_rxnA_0_F3_1.txt')
temp2 = np.loadtxt('chloro_ref_rxnC_0_F5_1.txt')

# normalize data set to the maximum intensity in the dataset
temp1[:,1] = temp1[:,1]/(np.amax(temp1[:,1]))
temp2[:,1] = temp2[:,1]/(np.amax(temp2[:,1]))

# determine local extrema; if you find that it picks too many peaks, increase
#   the value for "order"
max1 = argrelextrema(temp1[:,1], np.greater, order = 5000)
max2 = argrelextrema(temp2[:,1], np.greater, order = 5000)

#plot intensity versus m/z
plt.plot(temp1[:,0], temp1[:,1], color = "grey", linewidth = 1.0, label = 'H3K9R')
plt.plot(temp2[:,0], temp2[:,1], color = "k", linewidth = 1.0, label = 'H3K9RMe')

# enter the bounds to be used for picking the maxima to plot (m/z)
win1_max = 2300
win1_min = 2200

win2_max = 2400
win2_min = 2200

# pick extrema for each data set
for i in max1[0]:
    # Conditional statement to calculate maxima of points in a specfic range.
    # adds a m/z label for picked maxima
   if ((temp1[:,0][i] > win1_) and (temp1[:,0][i] < 2300 )):
       plt.annotate(str(temp1[:,0][i]), xy=(temp1[:,0][i], temp1[:,1][i]+0.1), horizontalalignment = 'center', fontsize = 8, rotation=45)


for i in max2[0]:
    # Conditional statement to calculate maxima of points in a specfic range.
    # adds a m/z label for picked maxima
    if ((temp2[:,0][i] > 2200 ) and (temp2[:,0][i] < 2400 )):
       plt.annotate(str(temp2[:,0][i]), xy=(temp2[:,0][i], temp2[:,1][i]+0.1), horizontalalignment = 'center', fontsize = 8, rotation=45)

## plotting stuff ##
plt.legend(loc='upper right')
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('H3K9R peptide methylation')
plt.xlim(2240, 2340)
plt.ylim([-0.05, 1.2])
plt.tight_layout()

# uncomment the following line to save spectrum plot
#plt.savefig('PRDM9_methylation.png', format = 'png', dpi = 300)
plt.show()
