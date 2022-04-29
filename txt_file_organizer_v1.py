#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 09:45:16 2021

Takes multiple .txt files input from fluorometer export and
re-saves the file sans header and footer // concatenates
all fluorescence intensity values against a single column
(index = 0) of wavelength values

@author: emeryusher
"""

import numpy as np

## I was overly ambitious here; you can use these input statements if you want
# sample_name = input('Enter the sample name prefix: ')
# urea_num = int(input('Enter the number of urea concentrations: '))
# data_start = int(input('Enter the line number where the data begins: '))
# data_end = int(input('Enter the line number for the last line of data: '))

# in the autosave program on the fluorometer, each file gets the same prefix
sample_name = '210614_W131G_urea'
# urea_num is the number of files/datapoints
urea_num = 16

# denotes the part of the text file that contains float values
# literally just the lines in the file where the data begin and end
data_start = 20
data_end = 70

# calculate the rows to skip and the max number of rows
sr = data_start - 1
mr = data_end - sr

# loop reads each text file and re-saves the file without the header and footer
for n in range(1, urea_num + 1):
    temp = np.loadtxt(str(sample_name)+ '-' +str(n)+ '.txt', skiprows = sr, max_rows = mr)

    np.savetxt(str(sample_name)+ '-' +str(n)+ '_parsed.txt', temp, '%1.3f')

temp_rest = []
x = 0

# this loop takes the parsed data files from previous loop and, for the first file,
#   begins the array with index 0 = wavelength, index 1 = intensity.
# then ONLY the intensity values from the remaining files are appended.
for i in range(1, urea_num + 1):
    temp = np.loadtxt(str(sample_name)+ '-' +str(i)+ '_parsed.txt')


    if i == 1:
        temp_one = np.array(temp)
        # np.savetxt('temp_all_test.txt', temp_one, '%1.3f')

    else:
        temp_rest.append(temp.T[1])

    x = x + 1

# array manipulation for my sanity
output = np.vstack(temp_rest)
melt_data = np.vstack((temp_one.T, temp_rest))
trans_melt = melt_data.T

# save to master formatted txt file with sample name prefix
np.savetxt(str(sample_name)+ '_formatted.txt', trans_melt, '%1.3f')
