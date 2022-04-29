#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 18:18:30 2021

This script takes input csv file(s) containing a table of data downloaded from
the Protein Data Bank (PDB), and parses the entries and counts the number of
each amino acid. Generates a bar plot of the frequency of each amino acid in the
folded protein dataset.

The output text file from this script ("pdb_AAs.txt") and the disordered protein
output data ("disprot_AAs.txt") are used by the script "idp_AA_enrinchment.py"
to calculate the fold enrichment of AAs in IDP sequences compared to folded
protein sequences.

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# array containing strings of file names for the csv data to be imported
file_array = ['rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_00001-02500.csv',
              'rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_02501-05000.csv',
              'rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_05001-07500.csv',
              'rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_07501-10000.csv',
              'rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_10001-12500.csv',
              'rcsb_pdb_custom_report_065d03b1e291bd0bf71484b5c30fcfc8_12501-13121.csv']


# use a loop to read and then append the data into a single array ("appended_df")
# (because the PDB only lets you download like 2000 entries at a time)
appended_df = []
for i in file_array:

    # the csv file headers are as follows: entry ID / PDB ID / Sequence /
    # Polymer Entity Sequence Length / Macromolecule Name / Accession Code(s)

    # use pandas to select desired columns and append only those in the array
    pdb_df = pd.read_csv(str(i), usecols = ['PDB ID', 'Sequence', 'Accession Code(s)'])

    appended_df.append(pdb_df)

master_df = pd.concat(appended_df)

# get the unique uniprot IDs from the list; this can be used later
acc_nos = np.array(master_df['Accession Code(s)'].unique())

# remove the lines that are duplicates, but preserve the first instance
#   this makes it so that each protein is only counted once, even if there are
#   multiple PDB deposits for a single molecule
clean_df = master_df.drop_duplicates(subset=['Accession Code(s)'])

# remove any sequences that contain 6X His tag in the structure (to prevent
#   overrepresentation of His)
clean_df2 = clean_df[clean_df["Sequence"].str.contains("HHHHHH") == False]

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

aa_strings = clean_df2['Sequence'].unique()

aa_lab = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# loop through each protein sequence and count the each amino acid type
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

# throw those numbers into an alphabetized array
aa = np.array([A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y])

# uncomment to save the AA counts to a text file:
#np.savetxt('pdb_AAs.txt', aa, fmt = '%s')

# determine the frequency of each AA in the folded dataset
total = sum(aa)
freq = aa / total

# save the accession numbers, if desired
#np.savetxt('uniprot_IDs.txt', acc_nos, fmt = '%s')

## plotting stuff ##
ypos = np.arange(len(aa_lab))

plt.bar(ypos, freq, color = 'grey', edgecolor = 'white', lw = 1)
plt.xticks(ypos, aa_lab)

plt.xlabel('amino acid')
plt.ylabel('frequency')

#plt.savefig('disprot_AAs.ps', format = 'ps', dpi = 300)

plt.show()
