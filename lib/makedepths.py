#! /usr/bin/env python

import sys, os
import numpy as np
import scipy.interpolate as si
import astropy.constants as const
sys.path.append('../../BART/code/')
import readtransit as rt
import kurucz_inten as ki
import wine


'''
Notes
planetradius --> planetrad
Wtemp --> temp
eclipseduration --> ecltime
filterfiles->filters
'''


def noiseup(kuruczfile, planetfile, filters, geometry, outdir, outpre, 
            mass=0.806, radius=0.756, temp=5000., planetrad=1.138, 
            distance=19.3, diameter=650., ecltime=1.827, seed=0):
    """
    This function takes an input star and planet spectrum and computes 
    transit/eclipse depths with photon noise, as observed by a 
    JWST-like telescope. 

    Inputs
    ------
    kuruczfile: string. Path/to/file of the Kurucz stellar model.
    planetfile: string. Path/to/file of the planet's spectrum.
    filters   : list, strings. Paths/to/filter files.
    geometry  : string. Viewing geometry. 'eclipse' or 'transit'
    outdir    : string. Path/to/directory/ where the outputs will be saved.
    outpre    : string. Prefix to use for all saved out files.
    mass      : float.  Mass   of the host star [M_sun]
    radius    : float.  Radius of the host star [R_sun]
    temp      : float.  Temperature of the host star [K]
    planetrad : float.  Radius of the planet [R_jup]
    distance  : float.  Distance to planetary system [pc]
    diameter  : float.  Telescope diameter [cm]
    ecltime   : float.  Duration of the secondary eclipse [hours]
    seed      : int.    Seed for random number generation.

    Outputs
    -------
    Four files are produced:
    - 
    - 
    - 
    - 

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
    distance  *= pc2cm
    ecltime   *= 3600 #hours --> seconds

    # Temperature and log(g) of star
    logg = np.log10(G * mass / radius**2) #log g (cm/s2)

    # Read and interpolate Kurucz grid, planet spectrum
    starfl, starwn, tmodel, gmodel = wine.readkurucz(kuruczfile, temp, logg)
    planetwn, planetfl             = rt.readspectrum(planetfile)

    nifilter  = []
    wnindices = []

    # Multiply by 4pi steradians, surface area
    # Planetary spectrum is already integrated over steradians
    sPower =   starfl * (4 * np.pi * radius**2)
    pPower = planetfl * (4 * np.pi * planetrad**2)

    # Multiply by eclipse duration
    sEnergy = sPower * ecltime
    pEnergy = pPower * ecltime

    # Spread out over distance
    sEdensity = sEnergy / (4 * np.pi * distance**2)
    pEdensity = pEnergy / (4 * np.pi * distance**2)

    # Multiply by telescope area (total energy received)
    sSED = sEdensity * np.pi * (diameter/2.)**2.
    pSED = pEdensity * np.pi * (diameter/2.)**2.

    # Interpolate stellar SED
    sSEDinterp = si.interp1d(starwn, sSED)
    isSED      = sSEDinterp(planetwn)

    if geometry == 'transit':
        # Interpolate stellar flux
        starflinterp = si.interp1d(starwn, starfl)
        istarfl      = starflinterp(planetwn)

    # Initialize arrays for band-integrated energy and mean wavelength
    bandintegratedstar   = np.zeros(len(filters))
    bandintegratedplanet = np.zeros(len(filters))
    meanwn               = np.zeros(len(filters))

    # Read filters. Resample to planetwn
    for i in np.arange(len(filters)):
        filtwn, filttransm   = wine.readfilter(filters[i])
        meanwn[i]            = np.mean(filtwn)
        nifilt, rsSED, wnind = wine.resample(planetwn, filtwn, filttransm, 
                                             starwn, sSED)
        nifilter.append(nifilt)
        wnindices.append(wnind)

    # Loop over each filter file, read the file, interpolate to the
    # wavelength array, and integrate over the bandpass, weighted by the
    # filter. 
    for i in np.arange(len(filters)):
        # integrate over the bandpass
        bandintegratedstar[i]       = wine.bandintegrate(
                                           isSED[wnindices[i][0]], planetwn, 
                                           nifilter[i], wnindices[i])
        if   geometry == 'eclipse':
            bandintegratedplanet[i] = wine.bandintegrate(
                                           pSED[wnindices[i][0]], 
                                           planetwn, nifilter[i], wnindices[i])
        elif geometry == 'transit':
            bandintegratedplanet[i] = wine.bandintegrate(
                                           planetfl[wnindices[i][0]], 
                                           planetwn, nifilter[i], wnindices[i])
        else:
            print("Invalid `geometry` specification.\n")
            sys.exit()

    # Divide by photon energy to get number of photons (counts)
    # Find total photon signal
    sphotons = bandintegratedstar   / (h * meanwn * c)
    if geometry == 'eclipse':
        pphotons = bandintegratedplanet / (h * meanwn * c)
        phot_tot = sphotons + pphotons
    else:
        phot_tot = sphotons #planet flux is negligible in transit geometry
        # Multiply by 1 minus the transit depth = signal during transit
        phot_tot_rat = phot_tot * (1. - bandintegratedplanet)

    # Noise it up
    poisson_tot = phot_tot**(.5)
    poisson_s   = sphotons**(.5)
    if geometry == 'eclipse':
        poisson_p = pphotons**(.5)
        noise     = np.random.normal(0, poisson_tot)
        # Add noise to the band-integrated photon counts
        noisedpts    = phot_tot + noise
        noisedplanet = pphotons + noise
        noisedstar   = sphotons + noise
        # Calculate eclipse depths
        depths = noisedplanet / phot_tot
        # Calculate eclipse depth uncertainty
        unc    = phot_tot/sphotons * ((poisson_tot/phot_tot)**2 + \
                                      (poisson_s  /sphotons)**2)**(.5)
    else:
        poisson_tot_rat = (phot_tot_rat)**(.5)
        # Propagate error
        unc    =      phot_tot_rat / phot_tot * \
                 ((poisson_tot_rat / phot_tot_rat)**2 + \
                  (poisson_tot     / phot_tot    )**2)**(.5)
        noise  = np.random.normal(0, poisson_tot_rat) / phot_tot_rat
        depths = bandintegratedplanet + noise

    # Save out files
    # Save filter file paths for BART
    with open(outdir+outpre+'filters.txt', 'w') as f:
        for i in range(len(filters)):
            f.write('../00inputs/filters/' + \
                    os.path.basename(filters[i]) + '\n')
    if geometry == 'eclipse':
        # Save depths to a file
        with open(outdir+outpre+'ecldepths.txt', 'w') as f:
            for i in range(len(depths)):
                f.write(str(depths[i]) + '\n')
        # Save uncertainties to a file
        with open(outdir+outpre+'ecluncs.txt', 'w') as f:
            for i in range(len(unc)):
                f.write(str(unc[i]) + '\n')
        # Save noiseless depths
        with open(outdir+outpre+'ecl_noiseless.txt', 'w') as f:
            for i in range(len(filters)):
                f.write(str(pphotons[i]/phot_tot[i]) + '\n')
    else:
        # Save depths to a file
        with open(outdir+outpre+'tradepths.txt', 'w') as f:
            for i in range(len(depths)):
                f.write(str(depths[i]) + '\n')
        # Save uncertainties to a file
        with open(outdir+outpre+'trauncs.txt', 'w') as f:
            for i in range(len(unc)):
                f.write(str(unc[i]) + '\n')
        # Save noiseless depths
        with open(outdir+outpre+'tra_noiseless.txt', 'w') as f:
            for i in range(len(filters)):
                f.write(str(bandintegratedplanet[i]) + '\n')


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
    noiseup(kuruczfile, iso_ecl,   filters, 'eclipse', s01, 'iso_')
    noiseup(kuruczfile, iso_tra,   filters, 'transit', s01, 'iso_')
    noiseup(kuruczfile, noinv_ecl, filters, 'eclipse', s02, 'noinv_')
    noiseup(kuruczfile, noinv_tra, filters, 'transit', s02, 'noinv_')
    noiseup(kuruczfile, inv_ecl,   filters, 'eclipse', s03, 'inv_')
    noiseup(kuruczfile, inv_tra,   filters, 'transit', s03, 'inv_')





