#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 12:55:18 2022

This script uses the output text files containing amino acid counts for disordered
proteins ("disprot_AAs.txt") and folded proteins ("pdb_AAs.txt") to determine
the AA enrichment in IDPs vs. folded proteins.

@author: emeryusher
"""

import numpy as np
import matplotlib.pyplot as plt

disordered = np.loadtxt('disprot_AAs.txt')
folded = np.loadtxt('pdb_AAs.txt')

aa_lab = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

disordered_tot = sum(disordered)
disordered_freq = disordered / disordered_tot

folded_tot = sum(folded)
folded_freq = folded / folded_tot

IDP_enrich = (disordered_freq / folded_freq) - 1

## plotting stuff ##
ypos = np.arange(len(aa_lab))

plt.bar(ypos, IDP_enrich, color = 'grey', edgecolor = 'white', lw = 1)

plt.xticks(ypos, aa_lab)

plt.xlabel('amino acid')
plt.ylabel('enrichment wrt folded proteins')
plt.savefig('idp_AA_enrichment.ps', format = 'ps', dpi = 300)

plt.show()
