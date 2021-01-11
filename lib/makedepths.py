#! /usr/bin/env python

import sys, os
import numpy as np
import scipy.interpolate as si
import astropy.constants as const
sys.path.append('../../BART/code/')
import readtransit as rt
import kurucz_inten as ki
import wine


def noiseup(snr, planetfile, filters, geometry, outdir, outpre, 
            kuruczfile, mass=0.806, radius=0.756, temp=5000., planetrad=1.138, 
            seed=0):
    """
    This function takes an input star and planet spectrum and computes 
    transit/eclipse depths with photon noise, as observed by a 
    JWST-like telescope. 

    Inputs
    ------
    snr       : float.  Desired SNR value.
    planetfile: string. Path/to/file of the planet's spectrum.
    filters   : list, strings. Paths/to/filter files.
    geometry  : string. Viewing geometry. 'eclipse' or 'transit'
    outdir    : string. Path/to/directory/ where the outputs will be saved.
    outpre    : string. Prefix to use for all saved out files.
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

    # Make sure `outdir` has trailing slash:
    if outdir[-1] != '/':
        outdir = outdir + '/'

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
    planetwn, planetfl             = rt.readspectrum(planetfile)



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

    # Uncertainty
    unc    = bandflux / snr

    # Save out files
    # Save filter file paths for BART
    with open(outdir+outpre+'filters.txt', 'w') as f:
        for i in range(len(filters)):
            f.write('../00inputs/filters/' + \
                    os.path.basename(filters[i]) + '\n')
    if   geometry == 'eclipse':
        strname = 'ecl'
    elif geometry == 'transit':
        strname = 'tra'
    # Save noised depths to a file
    with open(outdir+outpre+strname+'depths.txt', 'w') as f:
        for i in range(len(depths)):
            f.write(str(bandflux[i]) + '\n')
    # Save uncertainties to a file
    with open(outdir+outpre+strname+'uncs.txt', 'w') as f:
        for i in range(len(unc)):
            f.write(str(unc[i]) + '\n')


if __name__ == '__main__':
    # Set paths
    kuruczfile = '../tests/00inputs/hd189733b-fp00k2odfnew.pck'
    filters    = ['../tests/00inputs/filters/0' + str(i) + '.dat' 
                  for i in range(53, 100)]

    iso_ecl   = '../code-output/01BART/c01hjcleariso/'   + \
                'iso_emission_spectrum.dat'
    iso_tra   = '../code-output/01BART/c01hjcleariso/'   + \
                'iso_transmission_spectrum.dat'
    noinv_ecl = '../code-output/01BART/c02hjclearnoinv/' + \
                'noinv_emission_spectrum.dat'
    noinv_tra = '../code-output/01BART/c02hjclearnoinv/' + \
                'noinv_transmission_spectrum.dat'
    inv_ecl   = '../code-output/01BART/c03hjclearinv/'   + \
                'inv_emission_spectrum.dat'
    inv_tra   = '../code-output/01BART/c03hjclearinv/'   + \
                'inv_transmission_spectrum.dat'

    s01 = '../tests/s01hjcleariso/'
    s02 = '../tests/s02hjclearnoinv/'
    s03 = '../tests/s03hjclearinv/'

    # Call the above function for all 6 synthetic cases
    snr = 50.0
    #noiseup(snr, iso_ecl,   filters, 'eclipse', s01, 'iso_', kuruczfile)
    #noiseup(snr, iso_tra,   filters, 'transit', s01, 'iso_', kuruczfile)
    #noiseup(snr, noinv_ecl, filters, 'eclipse', s02, 'noinv_', kuruczfile)
    #noiseup(snr, noinv_tra, filters, 'transit', s02, 'noinv_', kuruczfile)
    #noiseup(snr, inv_ecl,   filters, 'eclipse', s03, 'inv_', kuruczfile)
    #noiseup(snr, inv_tra,   filters, 'transit', s03, 'inv_', kuruczfile)

    # Uniform abundance cases
    iso_uni_ecl   = '../code-output/01BART/c01hjcleariso/'   + \
                    'iso_uni_emission_spectrum.dat'
    iso_uni_tra   = '../code-output/01BART/c01hjcleariso/'   + \
                    'iso_uni_transmission_spectrum.dat'
    noinv_uni_ecl = '../code-output/01BART/c02hjclearnoinv/' + \
                    'noinv_uni_emission_spectrum.dat'
    noinv_uni_tra = '../code-output/01BART/c02hjclearnoinv/' + \
                    'noinv_uni_transmission_spectrum.dat'
    inv_uni_ecl   = '../code-output/01BART/c03hjclearinv/'   + \
                    'inv_uni_emission_spectrum.dat'
    inv_uni_tra   = '../code-output/01BART/c03hjclearinv/'   + \
                    'inv_uni_transmission_spectrum.dat'

    s01 = '../tests/s01hjcleariso/'
    s02 = '../tests/s02hjclearnoinv/'
    s03 = '../tests/s03hjclearinv/'

    # Call the above function for all 6 synthetic cases
    snr_ecl =  50.0
    snr_tra = 300.0
    noiseup(snr_ecl, iso_uni_ecl,   filters, 'eclipse', s01, 'iso_uni_', kuruczfile)
    noiseup(snr_tra, iso_uni_tra,   filters, 'transit', s01, 'iso_uni_', kuruczfile)
    noiseup(snr_ecl, noinv_uni_ecl, filters, 'eclipse', s02, 'noinv_uni_', kuruczfile)
    noiseup(snr_tra, noinv_uni_tra, filters, 'transit', s02, 'noinv_uni_', kuruczfile)
    noiseup(snr_ecl, inv_uni_ecl,   filters, 'eclipse', s03, 'inv_uni_', kuruczfile)
    noiseup(snr_tra, inv_uni_tra,   filters, 'transit', s03, 'inv_uni_', kuruczfile)


