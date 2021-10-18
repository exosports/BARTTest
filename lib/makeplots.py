#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append("../../BART/modules/transit/scripts/")
import readtransit as rt
import scipy.constants as const
import voigtcomp

"""
This file contains three functions used to produce plots of spectra, and 
a main function to produce plots for each test. An output directory is 
specified in main().

plotspectrum() is used to plot the entire spectrum produced by Transit.
plotspeczoom() is used to plot a small part of the spectrum.
plotspeciso()  is used to plot the isothermal spectrum with the Planck function 
               overplotted.
plotspecabun() is used to plot the varying-abundance spectra produced 
               in the abundance test.

Run this from BARTTest/

Plots will be saved into the directory where spectrum.dat file is located 
unless otherwise specified.
"""

def plotspectrum(fname, geo, wl=True, oname=False, title=False):
    """
    This function produces a plot of the spectrum produced by Transit.
    
    Inputs
    ------
    fname: string. File name of the spectrum file, with directory relative to 
                   /BARTTest/.
    geo  : string. Geometry of the produced spectrum.
                   'transit' or 'eclipse'
    wl   : bool.   True plots over wavelength, False plots over wavenumber.
    oname: string. If oname is False, the default plotting name is used when 
                   saving the file. If it is a string, then that 
                   path/to/plotname is used.
    title: bool.   True produces a plot with a title. False does not.
    
    Outputs
    -------
    PNG file of the spectrum. Default naming convention is i.e.
      oneline_emission_spectrum.png
      fewline_transmission_spectrum.png
    
    Example
    -------
    >>> plotspectrum('fewline/fewline_emission_spectrum.dat', 'eclipse')
    """
    # Load the data
    try:
        if wl==True:
            wlength, flux = rt.readspectrum(fname, 0)
        else:
            wlength, flux = rt.readspectrum(fname)
    except:
        return

    # Get test name
    testname = fname.split('/')[-1].split('_')[0]

    # Plot the data
    plt.figure(0, (8,5))
    plt.clf()
    plt.plot(wlength, flux, "b")
    
    # Set title and plot filename based on the geometry
    if geo=='eclipse':
        if title==True:
            plt.title(testname + " Emission Spectrum")
        plt.ylabel("Flux  (erg s$^{-1}$ cm$^{-1}$)")
        plotname = testname + '_emission_spectrum.png'
    else:
        if title==True:
            plt.title(testname + " Transmission Spectrum")
        plt.ylabel("Modulation  ($(R_p/R_s)^2$)")
        plotname = testname + '_transmission_spectrum.png'
    
    # Set X axis label to that specified
    if wl==True:
        plt.xlabel(u"Wavelength  (\u00b5m)")
    else:
        plt.xlabel("Wavenumber  (cm$^{-1}$)")

    # Set x limits
    plt.xlim(np.amin(wlength), np.amax(wlength))
    
    # Save the plot
    if oname==False:
        plt.savefig(fname.split('/')[0] + '/' + plotname)
    else:
        plt.savefig(oname)
    plt.clf()


def plotspeczoom(fname, geo, loc, wl=True, xlims=False, oname=False, 
                 title=True):
    """
    This function produces a plot of a portion of the spectrum produced by 
    Transit.
    
    Inputs
    ------
    fname: string. File name of the spectrum file.
    geo  : string. Geometry of the produced spectrum.
                   'transit' or 'eclipse'
    loc  : float.  Location in microns (if wl=True) or inverse centimeters 
                   (if wl=False) to zoom in on.
    wl   : bool.   True plots over wavelength, False plots over wavenumber.
    xlims: tuple.  (xmin, xmax) Minimum and maximum X-axis values for the plot. 
                   Only use this if user desires a different range for the 
                   X axis than the default.
    oname: string. If oname is False, the default plotting name is used. If 
                   it is a string, then that path/to/plotname is used.
    title: bool.   True plots w/ title. False does not.
        
    Outputs
    -------
    PNG file of the portion of the spectrum. Naming convention is i.e.
      fewline_emission_spectrum_zoom2500nm.png

    Example
    -------
    >>> plotspeczoom('fewline/fewline_emission_spectrum.dat', 'eclipse', 2.5)
    """
    # Load the data
    try:
        if wl==True:
            wlength, flux = rt.readspectrum(fname, 0)
            # File is ordered by decreasing wavelength--reverse it
            wlength = wlength[::-1]
            flux    = flux   [::-1]
        else:
            wlength, flux = rt.readspectrum(fname)
    except:
        return

    # Get test name
    testname = fname.split('/')[-1].split('_')[0]
    
    # Find where to zoom in, trim data to this region
    hi       = wlength[wlength < loc]
    wlentrim = wlength[len(hi) - 33 : len(hi) + 33] # This is why reverse wl
    fluxtrim = flux   [len(hi) - 33 : len(hi) + 33] # above
    
    # Plot the data
    plt.figure(0, (8,5))
    plt.clf()
    plt.plot(wlentrim, fluxtrim, "b")
    
    # Set title and plot filename based on the geometry
    if geo=='eclipse':
        if title==True:
            plt.title(testname + " Emission Spectrum")
        plt.ylabel("Flux  (erg s$^{-1}$ cm$^{-1}$)")
        plotname = testname + '_emission_spectrum_zoom'
    else:
        if title==True:
            plt.title(testname + " Transmission Spectrum")
        plt.ylabel("Modulation  ($(R_p/R_s)^2$)")
        plotname = testname + '_transmission_spectrum_zoom'
    
    # Set X axis label to that specified
    if wl==True:
        plt.xlabel(u"Wavelength  (\u00b5m)")
        plotname += str(1000*loc)[:4] + 'nm'
    else:
        plt.xlabel("Wavenumber  (cm$^{-1}$)")
        plotname += str(loc) + 'cm-1'
    
    # Set the X axis limits for the plot
    plt.xlim(wlentrim[0], wlentrim[-1])
    
    # Save the plot
    if oname==False:
        plt.savefig(fname.split('/')[0] + '/' + plotname + '.png')
    else:
        plt.savefig(oname)
    plt.clf()


def plotspeciso(fname, atm, geo, wl=True, oname=False, title=False):
    """
    This function produces a plot of the isothermal spectrum produced by 
    the isothermal test with the Planck function overplotted.
    
    Inputs
    ------
    fname: string. File name of the spectrum file, with directory relative to 
                   /BARTTest/.
    atm  : string. File name of the atmospheric file used, with directory.
    geo  : string. Geometry of the produced spectrum.
                   'transit' or 'eclipse'
    wl   : bool.   True plots over wavelength, False plots over wavenumber.
    oname: string. If oname is False, the default plotting name is used when 
                   saving the file. If it is a string, then that 
                   path/to/plotname is used.
    title: bool.   If True, plots w/ title. If False, it does not.
    
    Outputs
    -------
    PNG file of the spectrum. Default naming convention is i.e.
      isothermal_emission_spectrum.png
      isothermal_transmission_spectrum.png
    
    Example
    -------
    >>> plotspeciso('isothermal/isothermal_emission_spectrum.dat', \
                    'isothermal/isothermal.atm', 'eclipse')

    Notes
    -----
    For this plot to be produced correctly, an isothermal atmospheric file must 
    be supplied. There is a check for this in the function. A message will be 
    printed to the screen, and the code will continue to run. The resulting plot
    will however not contain the overplotted Planck function.

    """
    # Get the temperatures from the atm file
    foo = open(atm, 'r')
    lines = foo.readlines()
    lines = lines[13:] # trim the first 12 lines as they are headers

    # Array to hold temps
    temparr = np.zeros(len(lines), dtype=float)

    # Read in all the temps
    for i in range(len(lines)):
        temparr[i] = lines[i].split()[2]

    # Check that they are all the same--if not, yell at the user
    if len(np.unique(temparr)) != 1:
        print("The atmospheric file supplied is non-isothermal! " + \
              "Please supply an isothermal atmospheric file.\n")
        bgtemp = None
    else:
        bgtemp = np.unique(temparr)

    # Load the data
    try:
        if wl:
            wlength, flux = rt.readspectrum(fname, 0)
        else:
            wlength, flux = rt.readspectrum(fname)
    except:
        print("\nIsothermal plot: invalid specification for the data file name.")
        print("Problematic file:", fname, '\n')
        return

    # Get test name
    testname = fname.split('/')[-1].split('_')[0]

    # Plot the data
    plt.clf()
    fig1   = plt.figure(1)
    frame1 = fig1.add_axes((.1, .3, .8, .6))
    plt.plot(wlength, flux, color="k", lw=2, label='Transit')
    
    # Set title and plot filename based on the geometry
    if geo=='eclipse':
        if title:
            plt.title(testname + " Emission Spectrum")
        plt.ylabel("Flux  (erg s$^{-1}$ cm$^{-1}$)", fontsize=14)
        plotname = testname + '_emission_spectrum.png'
    else:
        if title:
            plt.title(testname + " Transmission Spectrum")
        plt.ylabel("Modulation  ($(R_p/R_s)^2$)", fontsize=14)
        plotname = testname + '_transmission_spectrum.png'
    
    # Plot the Planck function
    if wl:
        # Convert to wavenumber
        wavenum = 10000./wlength
        # Calculate the Planck function
        planck = 2. * const.h*1e7 * (wavenum**3) * (const.c*100)**2 /   \
                 (np.exp(const.h*1e7 * wavenum * const.c*100 /          \
                 (const.k*1e7) / bgtemp) - 1)
        # Plot it. Multiply Planck by pi because of how Transit calcs flux
        # Reverse Planck to match wlength
        plt.plot(wlength, np.pi*planck, color='goldenrod', lw=2, ls="--", 
                 label='Planck')
    else:
        planck = 2. * const.h*1e7 * (wlength**3) * (const.c*100)**2 /   \
                 (np.exp(const.h*1e7 * wlength * const.c*100 /          \
                 (const.k*1e7) / bgtemp) - 1)
        plt.plot(wlength, np.pi*planck, color='#ffff00', lw=2, ls="--", 
                 label='Planck')
    plt.legend(loc='upper right')

    frame1.set_xticklabels([])

    # Residuals
    frame2 = fig1.add_axes((.1, .1, .8, .2))
    plt.plot(wlength, 100*(flux - (np.pi*planck)) / (np.pi*planck))
    plt.ylabel('Residuals (%)', fontsize=14)
    yticks = frame2.yaxis.get_major_ticks()
    yticks[-1].label1.set_visible(False)

    # Set X axis label to that specified
    if wl:
        plt.xlabel(u"Wavelength  (\u00b5m)", fontsize=14)
    else:
        plt.xlabel("Wavenumber  (cm$^{-1}$)", fontsize=14)
    
    
    # Save the plot
    if not oname:
        plt.savefig(fname.split('/')[0] + '/' + plotname, bbox_inches='tight')
    else:
        plt.savefig(oname, bbox_inches='tight')
    plt.clf()



def plotspecabun(fnames, base, geo='eclipse', xlims=(2.28905,2.28945), 
                 oname=False, title=False, fext='.pdf'):
    """
    The function produces a plot of the spectra produced in abundance.
    All spectra are plotted together to show the difference in line depth.

    Inputs
    ------
    fnames: list of strings. File names of the spectrum files.
                   [fname1, fname2, fname3, ...]
    base  : string. File name of the spectrum with the line moved.
    geo   : string. Geometry of the produced spectrum.
                   'transit' or 'eclipse'
    xlims : tuple. Minimum and maximum X-axis values to be plotted.
                   Default corresponds to the default line used in abundance.
                   (xmin, xmax)
    oname: string. If oname is False, the default plotting name is used. If 
                   it is a string, then that path/to/plotname is used.
    title : bool   If True, plots w/ title. If False, it does not.

    Outputs
    -------
    PNG file of the portion of the spectrum. Naming convention is i.e.
      abundance_emission_spectra.png
    """
    # Load and plot spectrum w/ line moved
    wl, fl = rt.readspectrum(base, 0)
    plt.plot(wl, fl, label='No line', color=(0,0,0))

    # Loop over list of files
    for fn in range(len(fnames)):
        # Load data
        wlength, flux = rt.readspectrum(fnames[fn], 0)
        col = float(fn)/float(len(fnames))
        # Plot it, with label and no repeat colors
        plt.plot(wlength, flux, label=fnames[fn].split('/')[-1].split('_')[1], \
                 color=(col, 0, 1-col))

    # Set axis limits for the plot
    plt.xlim(xlims[0], xlims[1])
    loc = np.where(fl - flux == np.amax(fl - flux))[0][0]
    plt.ylim(flux[loc] - 5, fl[loc] + 5)
    
    # Label rest of plot, save
    if title:
        plt.title('Varying Abundance of One Line, \n Optically Thin Regime')
    plt.ylabel("Flux  (erg s$^{-1}$ cm$^{-1}$)", fontsize=14)
    plt.xlabel(u"Wavelength  (\u00b5m)", fontsize=14)
    plt.legend(loc="best")
    if not oname:
        plt.savefig(base.split('/')[0] + '/' + 'abundance_emission_spectra'+fext,
                    bbox_inches='tight')
    else:
        plt.savefig(oname+fext, bbox_inches='tight')
    plt.clf()


def main():
    """
    This function produces plots of the spectra produced in the tests.
    """
    print("Producing plots of spectra...")

    odir = '../code-output/01BART/'
    rdir = '../results/01BART/'
    
    # oneline
    try:
        plotspectrum(odir+'f01oneline/oneline_emission_spectrum.dat', 
                     'eclipse', 
                     oname=rdir+'f01oneline_emission_spectrum.png')
    except Exception as e:
        print(e)
        pass

    # fewline
    try:
        plotspectrum(odir+'f02fewline/fewline_emission_spectrum.dat', 
                     'eclipse', 
                     oname=rdir+'f02fewline_emission_spectrum.png')
        plotspectrum(odir+'f02fewline/fewline_transmission_spectrum.dat', 
                     'transit', 
                     oname=rdir+'f02fewline_transmission_spectrum.png')
    except Exception as e:
        print(e)
        pass

    # multiline
    try:
        plotspectrum(odir+'f03multiline/multiline_emission_spectrum.dat', 
                     'eclipse', 
                     oname=rdir+'f03multiline_emission_spectrum.png')
        plotspectrum(odir+'f03multiline/multiline_transmission_spectrum.dat', 
                     'transit', 
                     oname=rdir+'f03multiline_transmission_spectrum.png')
    except Exception as e:
        print(e)
        pass

    # broadening
    try:
        plotspectrum(odir+'f04broadening/broadening_emission_spectrum.dat', 
                     'eclipse', 
                     oname=rdir+'f04broadening_emission_spectrum.png')
        voigtcomp.comp()
    except Exception as e:
        print(e)
        pass

    # abundance
    fnames = [odir+'f05abundance/abundance_1e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_2e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_3e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_4e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_5e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_6e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_7e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_8e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_9e-4_emission_spectrum.dat',
              odir+'f05abundance/abundance_1e-3_emission_spectrum.dat']
    try:
        plotspecabun(fnames, 
                     odir+'f05abundance/abundance_0_emission_spectrum.dat', 
                     oname=rdir+'f05abundance_emission_spectra')
    except Exception as e:
        print(e)
        pass

    # blending
    try:
        plotspectrum(odir+'f06blending/blending_emission_spectrum.dat', 
                     'eclipse', 
                     oname=rdir+'f06blending_emission_spectrum.png')
    except Exception as e:
        print(e)
        pass

    # multicia
    try:
        plotspectrum(odir+'f07multicia/noCIA_emission_spectrum.dat',  'eclipse', 
                     oname=rdir+'f07noCIA_emission_spectrum.png')
        plotspectrum(odir+'f07multicia/oneCIA_emission_spectrum.dat', 'eclipse', 
                     oname=rdir+'f07oneCIA_emission_spectrum.png')
        plotspectrum(odir+'f07multicia/twoCIA_emission_spectrum.dat', 'eclipse', 
                     oname=rdir+'f07twoCIA_emission_spectrum.png')
    except Exception as e:
        print(e)
        pass
    
    # isothermal
    try:
        plotspeciso(odir+'f08isothermal/isothermal_emission_spectrum.dat', 
                    '../tests/f08isothermal/isothermal.atm', 'eclipse', 
                    oname=rdir+'f08isothermal_emission_spectrum.pdf')
    except Exception as e:
        print(e)
        pass

    


if __name__ == "__main__":
    main()

