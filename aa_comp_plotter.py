#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 18:18:30 2021

This script takes an excel (xlsx) file as an input and generates a plot of amino
acid usage frequency in IDPs. The dataset was obtained from the Disprot human
IDP database. Prior to importing here, the excel file was sorted such that the
longest IDR sequence was first (but this was mostly because I'm not good at pandas).

NB: any entries that contained sequences of less than 30 residues were omitted
from this analysis (per the convention for what length defines an IDP/IDR).

The output text file from this script ("disprot_AAs.txt") and the folded protein
output data ("pdb_AAs.txt") are used by the script "idp_AA_enrinchment.py"
to calculate the fold enrichment of AAs in IDP sequences compared to folded
protein sequences.

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# excel sheet was pre-sorted such that the entry with the longest IDR was first
#   NB: you can skip this step and sort via pandas instead if you are so inclined
disprot_df = pd.read_excel('disprot_Hsapiens_v2.xlsx', usecols = ['acc', 'disprot_id', 'region_sequence'])


# get the unique uniprot IDs from the list; this can be used later
acc_nos = np.array(disprot_df['acc'].unique())

# remove the lines that are duplicates, but preserve the first instance (which contains the longest IDR sequence)
clean_df = disprot_df.drop_duplicates(subset=['acc'])

# make variables for each amino acid that we will iteratively add to while counting
#   the appearnce of each residue over all of the entries

A = 0
C = 0
D = 0
E = 0
F = 0
G = 0
H = 0
I = 0
K = 0
L = 0
M = 0
N = 0
P = 0
Q = 0
R = 0
S = 0
T = 0
V = 0
W = 0
Y = 0

aa_strings = clean_df['region_sequence'].unique()

aa_lab = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

for i in range(0, len(aa_strings)):
    for x in aa_strings[i]:
        if x == 'A':
            A = A + 1
        elif x == 'C':
            C = C + 1
        elif x == 'D':
            D = D + 1
        elif x == 'E':
            E = E + 1
        elif x == 'F':
            F = F + 1
        elif x == 'G':
            G = G + 1
        elif x == 'H':
            H = H + 1
        elif x == 'I':
            I = I + 1
        elif x == 'K':
            K = K + 1
        elif x == 'L':
            L = L + 1
        elif x == 'M':
            M = M + 1
        elif x == 'N':
            N = N + 1
        elif x == 'P':
            P = P + 1
        elif x == 'Q':
            Q = Q + 1
        elif x == 'R':
            R = R + 1
        elif x == 'S':
            S = S + 1
        elif x == 'T':
            T = T + 1
        elif x == 'V':
            V = V + 1
        elif x == 'W':
            W = W + 1
        elif x == 'Y':
            Y = Y + 1

aa = np.array([A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y])

total = sum(aa)
freq = aa / total

# save amino acid counts in text file
np.savetxt('disprot_AAs.txt', aa, fmt = '%s')

## plotting stuff ##
ypos = np.arange(len(aa_lab))

plt.bar(ypos, freq, color = 'grey', edgecolor = 'white', lw = 1)
plt.xticks(ypos, aa_lab)

plt.xlabel('amino acid')
plt.ylabel('frequency')
#plt.savefig('disprot_AAs.ps', format = 'ps', dpi = 300)

plt.show()
