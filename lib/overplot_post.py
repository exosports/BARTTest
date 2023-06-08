#! /usr/bin/env python

import sys, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.ion()


def comp_histogram(stack1, stack2, name1, name2, 
                   parname=None, 
                   fignum=-12, fs=30, savefile=None, bins=60):
    """
    Plots the probability density functions of 1D marginalized posteriors 
    for two MCMC runs.
    
    Inputs
    ------
    stack1  : array.  Posterior, shaped (nparams, niterations)
    stack2  : array.  Same as `stack1`, but for the posterior to be compared.
    name1   : string. Label name for `stack1`
    name2   : string. Label name for `stack2`
    parname : list, strings. Parameter names.
    fs      : int.    Font size for plots.
    savefile: string. Path/to/file where the plot will be saved.
    bins    : int.    Number of bins for the histograms.
    """
    if np.shape(stack1)[0] != np.shape(stack2)[0]:
        raise ValueError('The posteriors must have the same ' + \
                         'number of parameters.')
    npars, niter1 = np.shape(stack1)
    npars, niter2 = np.shape(stack2)

    # Set default parameter names:
    if parname is None:
        namelen = int(2+np.log10(np.amax([npars-1,1])))
        parname = np.zeros(npars, "<U%d"%namelen)
        for i in np.arange(npars):
            parname[i] = "P" + str(i).zfill(namelen-1)

    # Set number of rows:
    if npars < 10:
        nperrow = 3
    else:
        nperrow = 4
    nrows = (npars - 1)//nperrow + 1
    # Set number of columns:
    if   npars > 9:
        ncolumns = 4
    elif npars > 4:
        ncolumns = 3
    else:
        ncolumns = (npars+2)//3 + (npars+2)%3  # (Trust me!)

    histheight = 4 + 4*(nrows)
    if nrows == 1:
        bottom = 0.25
    else:
        bottom = 0.15

    fig = plt.figure(fignum, figsize=(18, histheight))
    plt.clf()
    plt.subplots_adjust(left=0.1, right=0.95, bottom=bottom, top=0.9,
                        hspace=0.8, wspace=0.25)

    for i in np.arange(npars):
        ax = plt.subplot(nrows, ncolumns, i+1)
        a  = plt.xticks(size=fs-4, rotation=90)
        a  = plt.yticks(size=fs-4)
        if i%ncolumns == 0:
            plt.ylabel('Normalized PDF', size=fs)
        plt.xlabel(parname[i], size=fs)
        rng   = min(stack1[i].min(), stack2[i].min()), \
                max(stack1[i].max(), stack2[i].max())
        a = plt.hist(stack1[i], bins, range=rng, alpha=0.5, label=name1, density=True, 
                     color='b')
        a = plt.hist(stack2[i], bins, range=rng, alpha=0.5, label=name2, density=True, 
                     color='r')
        if i == npars - 1:
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={"size":fs})
    fig.align_labels()

    if savefile is not None:
        plt.savefig(savefile, bbox_inches='tight')

    plt.close()


if __name__ == '__main__':
    # Barstow et al (2020), NEMESIS
    grid10 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-nemesis-old/output.npy')
    grid01 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-nemesis-nocloud/output.npy')
    # Stack both
    stack10 = grid10[0]
    stack01 = grid01[0]
    for c in range(1, grid10.shape[0]):
        stack10 = np.hstack((stack10, grid10[c, :, 40000:])) #burnin of 40,000
    for c in range(1, grid01.shape[0]):
        stack01 = np.hstack((stack01, grid01[c, :, 50000:])) #burnin of 50,000
    # Shift the molecular abundances
    stack10[2:] -= 6
    stack01[2:] -= 6
    # Plot
    comp_histogram(stack01, stack10, 
                   '0.1 cm$^{-1}$ grid', 
                   '1.0 cm$^{-1}$ grid', 
                   parname=['T (K)',  'R$_p$ (km)',   'H$_2$O', 'CO'], 
                   savefile='../results/01BART/s04hjcleariso-BarstowEtal-model0-nemesis-postcompgrid.pdf')

    # Barstow et al (2020), CHIMERA
    grid10 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-chimera-old/output.npy')
    grid01 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-chimera-nocloud/output.npy')
    # Stack both
    stack10 = grid10[0]
    stack01 = grid01[0]
    for c in range(1, grid10.shape[0]):
        stack10 = np.hstack((stack10, grid10[c, :, 40000:])) #burnin of 40,000
    for c in range(1, grid01.shape[0]):
        stack01 = np.hstack((stack01, grid01[c, :, 50000:])) #burnin of 50,000
    # Shift the molecular abundances
    stack10[2:] -= 6
    stack01[2:] -= 6
    # Plot
    comp_histogram(stack01, stack10, 
                   '0.1 cm$^{-1}$ grid', 
                   '1.0 cm$^{-1}$ grid', 
                   parname=['T (K)',  'R$_p$ (km)',   'H$_2$O', 'CO'], 
                   savefile='../results/01BART/s04hjcleariso-BarstowEtal-model0-chimera-postcompgrid.pdf')

    # Barstow et al (2020), Tau-REx
    grid10 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-taurex-old/output.npy')
    grid01 = np.load('../code-output/01BART/s04hjcleariso-BarstowEtal-model0-taurex-nocloud/output.npy')
    # Stack both
    stack10 = grid10[0]
    stack01 = grid01[0]
    for c in range(1, grid10.shape[0]):
        stack10 = np.hstack((stack10, grid10[c, :, 40000:])) #burnin of 40,000
    for c in range(1, grid01.shape[0]):
        stack01 = np.hstack((stack01, grid01[c, :, 50000:])) #burnin of 50,000
    # Shift the molecular abundances
    stack10[2:] -= 6
    stack01[2:] -= 6
    # Plot
    comp_histogram(stack01, stack10, 
                   '0.1 cm$^{-1}$ grid', 
                   '1.0 cm$^{-1}$ grid', 
                   parname=['T (K)',  'R$_p$ (km)',   'H$_2$O', 'CO'], 
                   savefile='../results/01BART/s04hjcleariso-BarstowEtal-model0-taurex-postcompgrid.pdf')

    # HD 189733 b
    grid10 = np.load('../code-output/01BART/r01hd189733b_ref1bar_lessdata_alpha-fixed_15max_grid1.0/output.npy')
    grid01 = np.load('../code-output/01BART/r01hd189733b_ref1bar_lessdata_alpha-fixed_15max_grid0.1/output.npy')
    # Stack both
    stack10 = grid10[0]
    stack01 = grid01[0]
    for c in range(1, grid10.shape[0]):
        stack10 = np.hstack((stack10, grid10[c, :, 50000:])) #burnin of 50,000
    for c in range(1, grid01.shape[0]):
        stack01 = np.hstack((stack01, grid01[c, :, 50000:])) #burnin of 50,000
    # Shift the molecular abundances
    stack10[3:] -= 6
    stack01[3:] -= 6
    # Plot
    comp_histogram(stack01, stack10, 
                   'HD 189733 b data - 0.1 cm$^{-1}$ grid', 
                   'HD 189733 b data - 1.0 cm$^{-1}$ grid', 
                   parname=['$\kappa$', '$\gamma_1$', "$\\beta$", 'H$_2$O', 'CO$_2$', 'CO', 'CH$_4$'], 
                   savefile='../results/01BART/r01hd189733b-postcompgrid.pdf')


