#! /usr/bin/env python

"""
This script builds the filter files for the IRS data as analyzed by Grillmair et al. (2008)

The mean wavelengths were obtained from Mike Line (priv. comm.), who obtained them from Carl Grillmair

We assume that the filters are symmetric about the mean wavelengths
The 2nd order data (5 - 8 microns) is spaced roughly evenly by ~0.12097 microns
The 1st order data (> 7.5 microns) is spaced roughly evenly by ~0.24193 microns
Note: The last 2nd order data point (mean wavelength of ~7.44 microns) therefore has a slight overlap 
      in the bandpass of the first 1st order data point (mean wavelength of ~7.55 microns)
"""

import sys, os
import numpy as np

wlen = np.array([5.26261, 5.38358, 5.50454, 5.62551, 5.74648, 5.86745, 5.98842, 6.10939, 6.23035, 6.35132, 6.47229, 6.59326, 6.71422, 6.83519, 6.95616, 7.07713, 7.1981, 7.31907, 7.44003, 
                 7.54588, 7.78781, 8.02975, 8.27168, 8.51361, 8.75555, 8.99749, 9.23942, 9.48135, 9.72329, 9.96523, 10.2072, 10.4491, 10.691, 10.933, 11.1749, 11.4168, 11.6588, 11.9007, 12.1426, 12.3846, 12.6265, 12.8684, 13.1104, 13.3523, 13.5942, 13.8362, 14.0781])

eps = 1e-8 # Small value, to be used when defining filter bandpasses

ordercut = 7.5

i2nd = wlen <  ordercut
i1st = wlen >= ordercut

# Arrays to hold the left/right endpoints of each bandpass
left2nd = np.zeros(np.sum(i2nd))
rght2nd = np.zeros(np.sum(i2nd))

left1st = np.zeros(np.sum(i1st))
rght1st = np.zeros(np.sum(i1st))

# Handle edge cases
# Assume width is defined by the midpoint to the nearest wavelength point
left2nd[ 0] = wlen[ 0] - (wlen[ 1] - wlen[ 0])/2.
rght2nd[ 0] = wlen[ 0] + (wlen[ 1] - wlen[ 0])/2.
left2nd[-1] = wlen[18] - (wlen[18] - wlen[17])/2.
rght2nd[-1] = wlen[18] + (wlen[18] - wlen[17])/2.

left1st[ 0] = wlen[19] - (wlen[20] - wlen[19])/2.
rght1st[ 0] = wlen[19] + (wlen[20] - wlen[19])/2.
left1st[-1] = wlen[-1] - (wlen[-1] - wlen[-2])/2.
rght1st[-1] = wlen[-1] + (wlen[-1] - wlen[-2])/2.

# Fill in the rest
for i in range(1, np.sum(i2nd)-1):
    left2nd[i] = wlen[i] - (wlen[i+1] - wlen[i])/2.
    rght2nd[i] = wlen[i] + (wlen[i+1] - wlen[i])/2.

for i in range(1, np.sum(i1st)-1):
    left1st[i] = wlen[19+i] - (wlen[19+i+1] - wlen[19+i])/2.
    rght1st[i] = wlen[19+i] + (wlen[19+i+1] - wlen[19+i])/2.

# Build the filters
basename = 'IRS_GrillmairEtal2008_'
fext     = '.dat'
for i in range(np.sum(i2nd)):
    with open(basename+str(i+1).zfill(2)+fext, 'w') as foo:
        foo.write(str(np.round(left2nd[i]-eps*2, 9))           + ' 0.000\n')
        foo.write(str(np.round(left2nd[i]-eps, 9))             + ' 0.000\n')
        foo.write(str(np.round(left2nd[i], 9))                 + ' 1.000\n')
        foo.write(str(np.round((left2nd[i]+rght2nd[i])/2., 9)) + ' 1.000\n')
        foo.write(str(np.round(rght2nd[i], 9))                 + ' 1.000\n')
        foo.write(str(np.round(rght2nd[i]+eps, 9))             + ' 0.000\n')
        foo.write(str(np.round(rght2nd[i]+eps*2, 9))           + ' 0.000\n')

for i in range(np.sum(i1st)):
    with open(basename+str(19+i+1).zfill(2)+fext, 'w') as foo:
        foo.write(str(np.round(left1st[i]-eps*2, 9))           + ' 0.000\n')
        foo.write(str(np.round(left1st[i]-eps, 9))             + ' 0.000\n')
        foo.write(str(np.round(left1st[i], 9))                 + ' 1.000\n')
        foo.write(str(np.round((left1st[i]+rght1st[i])/2., 9)) + ' 1.000\n')
        foo.write(str(np.round(rght1st[i], 9))                 + ' 1.000\n')
        foo.write(str(np.round(rght1st[i]+eps, 9))             + ' 0.000\n')
        foo.write(str(np.round(rght1st[i]+eps*2, 9))           + ' 0.000\n')


