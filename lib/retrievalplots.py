#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
sys.path.append("../../BART/code/")
import reader as rd
import constants as c
import bestFit as bf
import makeatm as ma
import PT as pt


plt.ion()


def readatm(atmfile):
    # Open the atm file and read it
    atmfoo = open(atmfile, 'r')
    lines  = atmfoo.readlines()

    # Store the molecule info
    mols = lines[9].split()

    # Find where the layer data starts
    for i in range(len(lines)):
        line = lines[i].split()
        try:
            if line[0] == '#Radius':
                break
        except:
            continue

    # Trim atm file
    lines = lines[i + 1:]

    # Array to hold atm file info
    atminfo = np.zeros((len(lines), len(mols) + 3), dtype=float)
    # Note that +3 is to hold radius, pressure, and temperature

    # Read in data
    for i in range(len(lines)):
        # Split info for layer `i`
        line = lines[i].split()
        for j in range(3 + len(mols)): # rad, press, temp, mol1, mol2, ...
            # Read in each column
            atminfo[i, j] = float(line[j])

    # Make sure atmosphere is ordered from bottom to top
    if atminfo[0, 0] > atminfo[1, 0]:
        for i in range(atminfo.shape[1]):
            atminfo[:, i] = atminfo[:, i][::-1]

    return mols, atminfo


def retrievedPT(datadir, atmfile, tepfile, nmol, solution, 
                outname, outdir=None, T_int=100.):
    """
    Inputs
    ------
    datadir:  string. Path/to/directory containing BART-formatted output.
    atmfile:  string. Path/to/atmospheric model file.
    tepfile:  string. Path/to/Transiting ExoPlanet file.
    nmol:     int.    Number of molecules being fit by MCMC.
    solution: string. Geometry of the system. 'eclipse' or 'transit'. 
    outname:  string. File name of resulting plot.
    outdir:   string. Path/to/dir to save `outname`. If None, defaults 
                      to the results directory of BARTTest if `datadir` is 
                      within BARTTest.
    T_int:    float.  Internal planetary temperature. Default is 100 K.
    """
    # Set outdir if not specified
    if outdir is None:
        try:
            if datadir[-1] != '/':
                datadir = datadir + '/'
            outdir = 'results'.join(datadir.rsplit('code-output',            \
                                                   1)).rsplit('/', 2)[0] + '/'
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
        except:
            print("Data directory not located within BARTTest.")
            print("Please specify an output directory `outdir` and try again.")
            sys.exit(1)

    # Read g_surf and R_planet from TEP file
    grav, Rp = ma.get_g(tepfile)

    # Read star data from TEP file, and semi-major axis
    R_star, T_star, sma, gstar = bf.get_starData(tepfile)

    # Read atmfile
    mols, atminfo = readatm(atmfile)
    pressure = atminfo[:, 1]

    # Read MCMC output file
    MCfile = datadir + 'MCMC.log'
    bestP, uncer = bf.read_MCMC_out(MCfile)
    allParams = bestP
    # Get number of burned iterations
    foo   = open(MCfile, 'r')
    lines = foo.readlines()
    foo.close()
    line = [foop for foop in lines if " Burned in iterations per chain:" in foop]
    burnin = int(line[0].split()[-1])

    # Figure out number of parameters
    nparams   = len(allParams)
    nradfit   = int(solution == 'transit')
    nPTparams = nparams - nmol - nradfit
    PTparams  = allParams[:nPTparams]

    # Plot the best PT profile
    kappa, gamma1, gamma2, alpha, beta = PTparams
    best_T = pt.PT_line(pressure, kappa,  gamma1, gamma2, alpha, beta, 
                        R_star,   T_star, T_int,  sma,    grav*1e2, 'const')

    # Load MCMC data
    MCMCdata = datadir + 'output.npy'
    data     = np.load(MCMCdata)
    nchains, npars, niter = np.shape(data)

    # Make datacube from MCMC data
    data_stack = data[0,:,burnin:]
    for c in np.arange(1, nchains):
        data_stack = np.hstack((data_stack, data[c, :, burnin:]))

    # Datacube of PT profiles
    PTprofiles = np.zeros((np.shape(data_stack)[1], len(pressure)))

    curr_PTparams = PTparams

    for k in np.arange(0, np.shape(data_stack)[1]):
        j = 0
        for i in np.arange(len(PTparams)):
            curr_PTparams[i] = data_stack[j,k]
            j +=1
        kappa, gamma1, gamma2, alpha, beta = curr_PTparams
        PTprofiles[k] = pt.PT_line(pressure, kappa,  gamma1, gamma2, 
                                   alpha, beta, R_star, T_star, T_int, 
                                   sma, grav*1e2, 'const')


    # Get percentiles (for 1, 2-sigma boundaries):
    low1 = np.percentile(PTprofiles, 16.0, axis=0)
    hi1  = np.percentile(PTprofiles, 84.0, axis=0)
    low2 = np.percentile(PTprofiles,  2.5, axis=0)
    hi2  = np.percentile(PTprofiles, 97.5, axis=0)
    median = np.median(PTprofiles, axis=0)

    # Plot and save figure
    plt.figure(2)
    plt.clf()
    ax=plt.subplot(111)
    ax.fill_betweenx(pressure, low2, hi2, facecolor="#62B1FF", edgecolor="0.5")
    ax.fill_betweenx(pressure, low1, hi1, facecolor="#1873CC",
                                                           edgecolor="#1873CC")
    plt.semilogy(median, pressure, "-", lw=2, label='Median',color="k")
    plt.semilogy(best_T, pressure, "-", lw=2, label="Best fit", color="r")
    plt.semilogy(atminfo[:,2], pressure, "--", lw=2, label='Input', color='r')
    plt.ylim(pressure[0], pressure[-1])
    plt.legend(loc="best")
    plt.xlabel("Temperature  (K)", size=15)
    plt.ylabel("Pressure  (bar)",  size=15)
    plt.savefig(outdir + outname)
    plt.close()


def retrievedabun(datadir, outname, atminput, 
                  atmretrv='bestFit.atm', outdir=None):
    """
    Inputs
    ------
    datadir:  string. Path/to/directory containing BART-formatted output.
    outname:  string. File name of resulting plot. Do not include path, use 
                      `outdir` for that purpose.
    atminput:  string. Path/to/atmospheric model file used for input.
    atmretrv:  string. Name of atmospheric model file retrieved with respect 
                       to `datadir`. Default is 'bestFit.atm', the default 
                       name of BART's output.
    outdir:   string. Path/to/dir to save `outname`. If None, defaults 
                      to the results directory of BARTTest.
    """
    # Make sure datadir has a trailing /
    if datadir[-1] != '/':
        datadir = datadir + '/'
    # Set outdir if not specified. If it is, ensure trailing /
    if outdir is None:
        try:
            outdir = 'results'.join(datadir.rsplit('code-output',            \
                                                   1)).rsplit('/', 2)[0] + '/'
        except:
            print("Data directory not located within BARTTest.")
            print("Please specify an output directory `outdir` and try again.")
            sys.exit(1)
    elif outdir[-1] != '/':
        outdir = outdir + '/'
    # Ensure outdir exists
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    # Read the atm files
    molin, atmin = readatm(atminput)
    molrs, atmrs = readatm(datadir+atmretrv)

    # Only worry about molecules that are being fit
    ignore = {'H', 'C', 'N', 'O', 'He', 'H2', 'N2'}

    # Colors for plotting
    cols = ['deep pink', 'medium purple', 'indian red', 'teal', 'gray', \
            'green', 'blue', 'red', 'purple', 'lime', 'black', 'maroon', 'aqua']

    # Plot abundances
    for i in range(len(molrs)):
        if molrs[i] in ignore:
            continue
        try:
            ind = molin.index(molrs[i])
        except:
            continue
        plt.plot(atmin[:,ind+3], atmin[:,1], ls='--', lw=2, \
                 color=cols[i], label='Input '+molin[ind])
        plt.plot(atmrs[:,i  +3], atmrs[:,1], ls='-', lw=2, color=cols[i], \
                 label='Retrieved '+molrs[i])

    # Labels, legend, save plot
    plt.xlabel('Molar Mixing Fraction')
    plt.xscale('log')
    plt.ylabel('Pressure [bars]')
    plt.yscale('log')
    plt.gca().invert_yaxis()
    plt.legend(loc='lower left', prop={'size': 9})
    plt.savefig(outdir+outname)
    plt.close()


if __name__ == '__main__':
    # Common inputs
    datadirbase = '../code-output/01BART/'
    testdirbase = '../tests/'
    tepfile     = '../tests/00inputs/HD189733b.tep'
    outdir      = '../results/01BART/'
    # Make the plots
    print('Plotting the PT profiles...\n')
    try:
        print('  s01hjcleariso...\n')
        retrievedPT(datadirbase+'s01hjcleariso-ecl/', 
                    testdirbase+'c01hjcleariso/iso.tea', 
                    tepfile, 5, 'eclipse', 's01hjcleariso-ecl-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedPT(datadirbase+'s01hjcleariso-tra/', 
                    testdirbase+'c01hjcleariso/iso.tea', 
                    tepfile, 5, 'transit', 's01hjcleariso-tra-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    try:
        print('  s02hjclearnoinv...\n')
        retrievedPT(datadirbase+'s02hjclearnoinv-ecl/', 
                    testdirbase+'c02hjclearnoinv/noinv.tea', 
                    tepfile, 5, 'eclipse', 's02hjclearnoinv-ecl-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedPT(datadirbase+'s02hjclearnoinv-tra/', 
                    testdirbase+'c02hjclearnoinv/noinv.tea', 
                    tepfile, 5, 'transit', 's02hjclearnoinv-tra-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    try:
        print('  s03hjclearinv...\n')
        retrievedPT(datadirbase+'s03hjclearinv-ecl/', 
                    testdirbase+'c03hjclearinv/inv.tea', 
                    tepfile, 5, 'eclipse', 's03hjclearinv-ecl-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedPT(datadirbase+'s03hjclearinv-tra/', 
                    testdirbase+'c03hjclearinv/inv.tea', 
                    tepfile, 5, 'transit', 's03hjclearinv-tra-PT.png', outdir)
    except Exception as e:
        print(e)
        pass

    print('Plotting the vertical abundance profiles...\n')

    try:
        print('  s01hjcleariso...\n')
        retrievedabun(datadirbase+'s01hjcleariso-ecl/', 
                      's01hjcleariso-ecl-abun.png', 
                      testdirbase+'c01hjcleariso/iso.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedabun(datadirbase+'s01hjcleariso-tra/', 
                      's01hjcleariso-tra-abun.png', 
                      testdirbase+'c01hjcleariso/iso.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass

    try:
        print('  s02hjclearnoinv...\n')
        retrievedabun(datadirbase+'s02hjclearnoinv-ecl/', 
                      's02hjclearnoinv-ecl-abun.png', 
                      testdirbase+'c02hjclearnoinv/noinv.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedabun(datadirbase+'s02hjclearnoinv-tra/', 
                      's02hjclearnoinv-tra-abun.png', 
                      testdirbase+'c02hjclearnoinv/noinv.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass

    try:
        print('  s03hjclearinv...\n')
        retrievedabun(datadirbase+'s03hjclearinv-ecl/', 
                      's03hjclearinv-ecl-abun.png', 
                      testdirbase+'c03hjclearinv/inv.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass

    try:
        retrievedabun(datadirbase+'s03hjclearinv-tra/', 
                      's03hjclearinv-tra-abun.png', 
                      testdirbase+'c03hjclearinv/inv.tea', 
                      outdir=outdir)
    except Exception as e:
        print(e)
        pass


