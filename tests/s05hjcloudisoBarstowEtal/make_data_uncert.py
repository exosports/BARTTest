#! /usr/bin/env python
"""
This is much like the makefilters*.py files in the s04 test
Since this test uses the same filters, this script only loads 
the Barstow et al. data and saves it out into files for convenience of 
copy/pasting into the BART configuration file.

Note: this script will not run unless you copy the relevant files from the 
      Barstow et al. compendium into this directory.
This script is meant as a reference for how the different BART inputs were 
constructed at the time of the test execution.
"""

import sys, os
import numpy as np


nemesis = np.loadtxt('nemesis_hd189733b-clouds_jwst_60.dat')
chimera = np.loadtxt('chimera_hd189733b-clouds_jwst_60.dat')
taurex  = np.loadtxt('taurex_hd189733b-clouds_jwst_60.dat')

np.savetxt('nemesis_data.txt', nemesis[:,1])
np.savetxt('nemesis_uncert.txt', nemesis[:,2])

np.savetxt('chimera_data.txt', chimera[:,1])
np.savetxt('chimera_uncert.txt', chimera[:,2])

np.savetxt('taurex_data.txt', taurex[:,1])
np.savetxt('taurex_uncert.txt', taurex[:,2])
