#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal      as ss
from scipy.ndimage.filters import gaussian_filter1d as gaussf
import astropy.constants as const

sys.path.append('../../BART/code/')
import readtransit as rt
import kurucz_inten as ki
import wine


plt.ion()

title  = False
fext   = '.pdf'
outdir = '../results/01BART/speccomp/'

kuruczfile = '../tests/00inputs/hd189733b-fp00k2odfnew.pck'

# Load the data
dat020 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.025grid_10max_spectrum.dat')
dat021 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.025grid_11max_spectrum.dat')
dat010 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.1grid_10max_spectrum.dat')
dat011 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.1grid_11max_spectrum.dat')
dat1 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.25grid_10max_spectrum.dat')
dat2 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_0.25grid_11max_spectrum.dat')
dat3 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_1.0grid_10max_spectrum.dat')
dat4 = np.loadtxt('../code-output/01BART/c02hjclearnoinv/noinv_uni_emission_1.0grid_11max_spectrum.dat')

# Load the filters
filters = ['../tests/00inputs/filters/0' + str(i) + '.dat' 
           for i in range(53, 100)]

# Unbinned, unsmoothed data over retrieval range
plt.plot(dat020[:,0], dat020[:,1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], dat021[:,1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], dat010[:,1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], dat011[:,1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], dat1[:,1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], dat2[:,1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], dat3[:,1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], dat4[:,1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Unsmoothed, unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_nosmooth'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], dat020[:,1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], dat021[:,1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], dat010[:,1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], dat011[:,1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Unsmoothed, unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_nosmooth_4small'+fext, bbox_inches='tight')
plt.close()

# Unbinned, smoothed data - Gaussian smoothing
# gaussf(spec, stdev)
g020_2  = gaussf(dat020[:,1],  2*40)
g020_5  = gaussf(dat020[:,1],  5*40)
g020_25 = gaussf(dat020[:,1], 25*40)

g021_2  = gaussf(dat021[:,1],  2*40)
g021_5  = gaussf(dat021[:,1],  5*40)
g021_25 = gaussf(dat021[:,1], 25*40)

g010_2  = gaussf(dat010[:,1],  2*10)
g010_5  = gaussf(dat010[:,1],  5*10)
g010_25 = gaussf(dat010[:,1], 25*10)

g011_2  = gaussf(dat011[:,1],  2*10)
g011_5  = gaussf(dat011[:,1],  5*10)
g011_25 = gaussf(dat011[:,1], 25*10)

g1_2  = gaussf(dat1[:,1],  2*4)
g1_5  = gaussf(dat1[:,1],  5*4)
g1_25 = gaussf(dat1[:,1], 25*4)

g2_2  = gaussf(dat2[:,1],  2*4)
g2_5  = gaussf(dat2[:,1],  5*4)
g2_25 = gaussf(dat2[:,1], 25*4)

g3_2  = gaussf(dat3[:,1],  2)
g3_5  = gaussf(dat3[:,1],  5)
g3_25 = gaussf(dat3[:,1], 25)

g4_2  = gaussf(dat4[:,1],  2)
g4_5  = gaussf(dat4[:,1],  5)
g4_25 = gaussf(dat4[:,1], 25)

plt.plot(dat020[:,0], g020_2, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_2, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_2, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_2, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], g1_2, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], g2_2, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], g3_2, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], g4_2, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=2 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-2'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], g020_5, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_5, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_5, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_5, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], g1_5, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], g2_5, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], g3_5, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], g4_5, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=5 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-5'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], g020_25, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_25, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_25, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_25, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], g1_25, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], g2_25, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], g3_25, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], g4_25, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=25 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-25'+fext, bbox_inches='tight')
plt.close()

# 4 smallest grid spacings
plt.plot(dat020[:,0], g020_2, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_2, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_2, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_2, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=2 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-2_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], g020_5, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_5, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_5, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_5, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=5 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-5_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], g020_25, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], g021_25, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], g010_25, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], g011_25, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=25 cm$^{-1}$), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_gsmooth-25_4small'+fext, bbox_inches='tight')
plt.close()

# Unbinned, smoothed data - Savitsky-Golay smoothing
# ss.savgol_filter(spec, nchannels, order)
sg020_5  = ss.savgol_filter(dat020[:,1],  5*40+1, 3)
sg020_25 = ss.savgol_filter(dat020[:,1], 25*40+1, 3)
sg020_75 = ss.savgol_filter(dat020[:,1], 75*40+1, 3)

sg021_5  = ss.savgol_filter(dat021[:,1],  5*40+1, 3)
sg021_25 = ss.savgol_filter(dat021[:,1], 25*40+1, 3)
sg021_75 = ss.savgol_filter(dat021[:,1], 75*40+1, 3)

sg010_5  = ss.savgol_filter(dat010[:,1],  5*10+1, 3)
sg010_25 = ss.savgol_filter(dat010[:,1], 25*10+1, 3)
sg010_75 = ss.savgol_filter(dat010[:,1], 75*10+1, 3)

sg011_5  = ss.savgol_filter(dat011[:,1],  5*10+1, 3)
sg011_25 = ss.savgol_filter(dat011[:,1], 25*10+1, 3)
sg011_75 = ss.savgol_filter(dat011[:,1], 75*10+1, 3)

sg1_5  = ss.savgol_filter(dat1[:,1],  5*4+1, 3)
sg1_25 = ss.savgol_filter(dat1[:,1], 25*4+1, 3)
sg1_75 = ss.savgol_filter(dat1[:,1], 75*4+1, 3)

sg2_5  = ss.savgol_filter(dat2[:,1],  5*4+1, 3)
sg2_25 = ss.savgol_filter(dat2[:,1], 25*4+1, 3)
sg2_75 = ss.savgol_filter(dat2[:,1], 75*4+1, 3)

sg3_5  = ss.savgol_filter(dat3[:,1],  5, 3)
sg3_25 = ss.savgol_filter(dat3[:,1], 25, 3)
sg3_75 = ss.savgol_filter(dat3[:,1], 75, 3)

sg4_5  = ss.savgol_filter(dat4[:,1],  5, 3)
sg4_25 = ss.savgol_filter(dat4[:,1], 25, 3)
sg4_75 = ss.savgol_filter(dat4[:,1], 75, 3)

plt.plot(dat020[:,0], sg020_5, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_5, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_5, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_5, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], sg1_5, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], sg2_5, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], sg3_5, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], sg4_5, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=5 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-5'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], sg020_25, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_25, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_25, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_25, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], sg1_25, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], sg2_25, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], sg3_25, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], sg4_25, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=25 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-25'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], sg020_75, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_75, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_75, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_75, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(dat1[:,0], sg1_75, label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(dat2[:,0], sg2_75, label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(dat3[:,0], sg3_75, label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(dat4[:,0], sg4_75, label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=75 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-75'+fext, bbox_inches='tight')
plt.close()

# 4 smallest bin spacings
plt.plot(dat020[:,0], sg020_5, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_5, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_5, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_5, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=5 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-5_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], sg020_25, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_25, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_25, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_25, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=25 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-25_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(dat020[:,0], sg020_75, label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(dat021[:,0], sg021_75, label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(dat010[:,0], sg010_75, label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(dat011[:,0], sg011_75, label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=75 cm$^{-1}$, cubic), unbinned')
plt.xlabel('Wavelength (um)')
plt.ylabel('Flux (erg/s/cm)')
plt.savefig(outdir+'spec_comp_nobin_savgolsmooth-75_4small'+fext, bbox_inches='tight')
plt.close()


### BINNED ###
def binit(planetwn, planetfl, filters, geometry, kuruczfile, 
          mass=0.806, radius=0.756, temp=5000., planetrad=1.138, 
          seed=0):
    """
    This function takes an input star and planet spectrum and computes 
    transit/eclipse depths with photon noise, as observed by a 
    JWST-like telescope. 

    Inputs
    ------
    planetfile: string. Path/to/file of the planet's spectrum.
    filters   : list, strings. Paths/to/filter files.
    geometry  : string. Viewing geometry. 'eclipse' or 'transit'
    kuruczfile: string. Path/to/file of the Kurucz stellar model.
    mass      : float.  Mass   of the host star [M_sun]
    radius    : float.  Radius of the host star [R_sun]
    temp      : float.  Temperature of the host star [K]
    planetrad : float.  Radius of the planet [R_jup]
    seed      : int.    Seed for random number generation.

    Outputs
    -------
    Four text files are produced:
    - list of filers
    - noiseless spectrum
    - noised spectrum
    - uncertainties on noised spectrum

    Notes
    -----
    The default values are based on the HD 189733 system.
    The default `ecltime` value comes from 
    https://academic.oup.com/mnras/article/395/1/335/1746726/Secondary-radio-eclipse-of-the-transiting-planet

    Revisions
    ---------
    2017-11-15  Ryan        Original implementation
    2018-04-18  Michael     Added transit calculations
    2019-06-06  Michael     Reworked into function, incorporated into BARTTest
    2020-08-08  Michael     Simplified computation
    """
    # Set the random seed for reproducibility
    np.random.seed(seed=seed)

    # Constants in CGS units
    G     = const.G.cgs.value
    h     = const.h.cgs.value
    c     = const.c.cgs.value
    Msun  = const.M_sun.cgs.value
    Rsun  = const.R_sun.cgs.value
    pc2cm = const.pc.cgs.value
    Rjup  = const.R_jup.cgs.value

    # Convert system parameters
    mass      *= Msun
    radius    *= Rsun
    planetrad *= Rjup

    rprs = planetrad / radius #ratio of planet/star radii

    # Temperature and log(g) of star
    logg = np.log10(G * mass / radius**2) #log g (cm/s2)

    # Read and interpolate Kurucz grid, planet spectrum
    starfl, starwn, tmodel, gmodel = wine.readkurucz(kuruczfile, temp, logg)

    # Initialize arrays for band-integrated spectrum and mean wavelength
    bandflux  = np.zeros(len(filters))
    meanwn    = np.zeros(len(filters))
    nifilter  = []
    wnindices = []
    istarfl   = [] # interpolated stellar flux
    # Read filters. Resample to planetwn
    for i in np.arange(len(filters)):
        filtwn, filttransm   = wine.readfilter(filters[i])
        meanwn[i]            = np.mean(filtwn)
        nifilt, strfl, wnind = wine.resample(planetwn, filtwn, filttransm,
                                             starwn, starfl)
        nifilter.append(nifilt)
        wnindices.append(wnind)
        istarfl.append(strfl)

    # Loop over each filter file, read the file, interpolate to the
    # wavelength array, and integrate over the bandpass, weighted by the
    # filter. 
    for i in np.arange(len(filters)):
        # integrate over the bandpass
      if   geometry == "eclipse":
        fluxrat = (planetfl[wnindices[i]]/istarfl[i]) * rprs*rprs
        bandflux[i] = wine.bandintegrate(fluxrat, planetwn,
                                      nifilter[i], wnindices[i])
      elif geometry == "transit":
        bandflux[i] = wine.bandintegrate(planetfl[wnindices[i]], planetwn,
                                      nifilter[i], wnindices[i])

    return meanwn, bandflux


# Binned, unsmoothed data
b020  = binit(10000./dat020[:,0], dat020[:,1], filters, 'eclipse', kuruczfile)
b021  = binit(10000./dat021[:,0], dat021[:,1], filters, 'eclipse', kuruczfile)
b010  = binit(10000./dat010[:,0], dat010[:,1], filters, 'eclipse', kuruczfile)
b011  = binit(10000./dat011[:,0], dat011[:,1], filters, 'eclipse', kuruczfile)
b1    = binit(10000./dat1[:,0], dat1[:,1], filters, 'eclipse', kuruczfile)
b2    = binit(10000./dat2[:,0], dat2[:,1], filters, 'eclipse', kuruczfile)
b3    = binit(10000./dat3[:,0], dat3[:,1], filters, 'eclipse', kuruczfile)
b4    = binit(10000./dat4[:,0], dat4[:,1], filters, 'eclipse', kuruczfile)

np.save('../tests/s02hjclearnoinv/noinv_uni_hires_ecl_noiseless.npy', b010[1])
np.save('../tests/s02hjclearnoinv/noinv_uni_hires_ecluncs_50snr.npy', b010[1] / 50.)

plt.plot(10000./b020[0], b020[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021[0], b021[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010[0], b010[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011[0], b011[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1[0], b1[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2[0], b2[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3[0], b3[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4[0], b4[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Unsmoothed, binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_nosmooth'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020[0], b020[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021[0], b021[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010[0], b010[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011[0], b011[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Unsmoothed, binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_nosmooth_4small'+fext, bbox_inches='tight')
plt.close()


# Binned, smoothed - Gaussian
b020_g2  = binit(10000./dat020[:,0], g020_2, filters, 'eclipse', kuruczfile)
b020_g5  = binit(10000./dat020[:,0], g020_5, filters, 'eclipse', kuruczfile)
b020_g25 = binit(10000./dat020[:,0], g020_25, filters, 'eclipse', kuruczfile)

b021_g2  = binit(10000./dat021[:,0], g021_2, filters, 'eclipse', kuruczfile)
b021_g5  = binit(10000./dat021[:,0], g021_5, filters, 'eclipse', kuruczfile)
b021_g25 = binit(10000./dat021[:,0], g021_25, filters, 'eclipse', kuruczfile)

b010_g2  = binit(10000./dat010[:,0], g010_2, filters, 'eclipse', kuruczfile)
b010_g5  = binit(10000./dat010[:,0], g010_5, filters, 'eclipse', kuruczfile)
b010_g25 = binit(10000./dat010[:,0], g010_25, filters, 'eclipse', kuruczfile)

b011_g2  = binit(10000./dat011[:,0], g011_2, filters, 'eclipse', kuruczfile)
b011_g5  = binit(10000./dat011[:,0], g011_5, filters, 'eclipse', kuruczfile)
b011_g25 = binit(10000./dat011[:,0], g011_25, filters, 'eclipse', kuruczfile)

b1_g2  = binit(10000./dat1[:,0], g1_2, filters, 'eclipse', kuruczfile)
b1_g5  = binit(10000./dat1[:,0], g1_5, filters, 'eclipse', kuruczfile)
b1_g25 = binit(10000./dat1[:,0], g1_25, filters, 'eclipse', kuruczfile)

b2_g2  = binit(10000./dat2[:,0], g2_2, filters, 'eclipse', kuruczfile)
b2_g5  = binit(10000./dat2[:,0], g2_5, filters, 'eclipse', kuruczfile)
b2_g25 = binit(10000./dat2[:,0], g2_25, filters, 'eclipse', kuruczfile)

b3_g2  = binit(10000./dat3[:,0], g3_2, filters, 'eclipse', kuruczfile)
b3_g5  = binit(10000./dat3[:,0], g3_5, filters, 'eclipse', kuruczfile)
b3_g25 = binit(10000./dat3[:,0], g3_25, filters, 'eclipse', kuruczfile)

b4_g2  = binit(10000./dat4[:,0], g4_2, filters, 'eclipse', kuruczfile)
b4_g5  = binit(10000./dat4[:,0], g4_5, filters, 'eclipse', kuruczfile)
b4_g25 = binit(10000./dat4[:,0], g4_25, filters, 'eclipse', kuruczfile)


plt.plot(10000./b020_g2[0], b020_g2[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g2[0], b021_g2[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g2[0], b010_g2[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g2[0], b011_g2[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_g2[0], b1_g2[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_g2[0], b2_g2[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_g2[0], b3_g2[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_g2[0], b4_g2[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=2 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-2'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_g5[0], b020_g5[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g5[0], b021_g5[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g5[0], b010_g5[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g5[0], b011_g5[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_g5[0], b1_g5[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_g5[0], b2_g5[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_g5[0], b3_g5[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_g5[0], b4_g5[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=5 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-5'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_g25[0], b020_g25[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g25[0], b021_g25[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g25[0], b010_g25[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g25[0], b011_g25[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_g25[0], b1_g25[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_g25[0], b2_g25[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_g25[0], b3_g25[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_g25[0], b4_g25[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=25 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-25'+fext, bbox_inches='tight')
plt.close()

# 4 smallest bin spacings
plt.plot(10000./b020_g2[0], b020_g2[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g2[0], b021_g2[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g2[0], b010_g2[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g2[0], b011_g2[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=2 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-2_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_g5[0], b020_g5[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g5[0], b021_g5[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g5[0], b010_g5[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g5[0], b011_g5[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=5 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-5_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_g25[0], b020_g25[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_g25[0], b021_g25[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_g25[0], b010_g25[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_g25[0], b011_g25[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('Gaussian smoothing (std=25 cm$^{-1}$), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_gsmooth-25_4small'+fext, bbox_inches='tight')
plt.close()

# Binned, smoothed - Savitsky-Golay
b020_sg5  = binit(10000./dat020[:,0], sg020_5, filters, 'eclipse', kuruczfile)
b020_sg25 = binit(10000./dat020[:,0], sg020_25, filters, 'eclipse', kuruczfile)
b020_sg75 = binit(10000./dat020[:,0], sg020_75, filters, 'eclipse', kuruczfile)

b021_sg5  = binit(10000./dat021[:,0], sg021_5, filters, 'eclipse', kuruczfile)
b021_sg25 = binit(10000./dat021[:,0], sg021_25, filters, 'eclipse', kuruczfile)
b021_sg75 = binit(10000./dat021[:,0], sg021_75, filters, 'eclipse', kuruczfile)

b010_sg5  = binit(10000./dat010[:,0], sg010_5, filters, 'eclipse', kuruczfile)
b010_sg25 = binit(10000./dat010[:,0], sg010_25, filters, 'eclipse', kuruczfile)
b010_sg75 = binit(10000./dat010[:,0], sg010_75, filters, 'eclipse', kuruczfile)

b011_sg5  = binit(10000./dat011[:,0], sg011_5, filters, 'eclipse', kuruczfile)
b011_sg25 = binit(10000./dat011[:,0], sg011_25, filters, 'eclipse', kuruczfile)
b011_sg75 = binit(10000./dat011[:,0], sg011_75, filters, 'eclipse', kuruczfile)

b1_sg5  = binit(10000./dat1[:,0], sg1_5, filters, 'eclipse', kuruczfile)
b1_sg25 = binit(10000./dat1[:,0], sg1_25, filters, 'eclipse', kuruczfile)
b1_sg75 = binit(10000./dat1[:,0], sg1_75, filters, 'eclipse', kuruczfile)

b2_sg5  = binit(10000./dat2[:,0], sg2_5, filters, 'eclipse', kuruczfile)
b2_sg25 = binit(10000./dat2[:,0], sg2_25, filters, 'eclipse', kuruczfile)
b2_sg75 = binit(10000./dat2[:,0], sg2_75, filters, 'eclipse', kuruczfile)

b3_sg5  = binit(10000./dat3[:,0], sg3_5, filters, 'eclipse', kuruczfile)
b3_sg25 = binit(10000./dat3[:,0], sg3_25, filters, 'eclipse', kuruczfile)
b3_sg75 = binit(10000./dat3[:,0], sg3_75, filters, 'eclipse', kuruczfile)

b4_sg5  = binit(10000./dat4[:,0], sg4_5, filters, 'eclipse', kuruczfile)
b4_sg25 = binit(10000./dat4[:,0], sg4_25, filters, 'eclipse', kuruczfile)
b4_sg75 = binit(10000./dat4[:,0], sg4_75, filters, 'eclipse', kuruczfile)

plt.plot(10000./b020_sg5[0], b020_sg5[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg5[0], b021_sg5[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg5[0], b010_sg5[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg5[0], b011_sg5[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_sg5[0], b1_sg5[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_sg5[0], b2_sg5[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_sg5[0], b3_sg5[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_sg5[0], b4_sg5[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=5 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-5'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_sg25[0], b020_sg25[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg25[0], b021_sg25[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg25[0], b010_sg25[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg25[0], b011_sg25[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_sg25[0], b1_sg25[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_sg25[0], b2_sg25[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_sg25[0], b3_sg25[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_sg25[0], b4_sg25[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=25 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-25'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_sg75[0], b020_sg75[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg75[0], b021_sg75[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg75[0], b010_sg75[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg75[0], b011_sg75[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b1_sg75[0], b1_sg75[1], label='0.25 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b2_sg75[0], b2_sg75[1], label='0.25 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b3_sg75[0], b3_sg75[1], label='1.0 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b4_sg75[0], b4_sg75[1], label='1.0 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=75 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-75'+fext, bbox_inches='tight')
plt.close()

# 4 smallest bin spacings
plt.plot(10000./b020_sg5[0], b020_sg5[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg5[0], b021_sg5[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg5[0], b010_sg5[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg5[0], b011_sg5[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=5 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-5_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_sg25[0], b020_sg25[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg25[0], b021_sg25[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg25[0], b010_sg25[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg25[0], b011_sg25[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=25 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-25_4small'+fext, bbox_inches='tight')
plt.close()

plt.plot(10000./b020_sg75[0], b020_sg75[1], label='0.025 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b021_sg75[0], b021_sg75[1], label='0.025 cm$^{-1}$ grid, 11 um max')
plt.plot(10000./b010_sg75[0], b010_sg75[1], label='0.1 cm$^{-1}$ grid, 10 um max')
plt.plot(10000./b011_sg75[0], b011_sg75[1], label='0.1 cm$^{-1}$ grid, 11 um max')
plt.xlim(10000./5500., 10000./2000.)
plt.legend(loc='center left', prop={'size' : 8}, bbox_to_anchor=(1, 0.5))
if title:
    plt.title('SavGol smoothing (len=75 cm$^{-1}$, cubic), binned')
plt.xlabel('Wavelength (um)')
plt.ylabel('F$_p$/F$_s$')
plt.savefig(outdir+'spec_comp_bin_savgolsmooth-75_4small'+fext, bbox_inches='tight')
plt.close()


