import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
#----------------------------------------------------------------------
def ReadSpectrum(fname):
    f = h5py.File(fname, "r")
    epe = np.array(f['epe'])
    spect = np.array(f['spect'])
    f.close()
    return epe, spect
#----------------------------------------------------------------------
def PlotSpectrumPath(ebind,spect,fac=1.0,kps=[],kpoints=[],klabel=[]):
    Ry = 27.211386

    fig, ax = plt.subplots(figsize=(5,5))

    kmin = 0
    kmax = 1
    if len(kps) > 0:
        kmin = np.amin(kps)
        kmax = np.amax(kps)


    smax = fac * np.amax(spect)
    wmin, wmax = Ry * np.amin(ebind), Ry * np.amax(ebind)
    asp =  (kmax - kmin) / (wmax - wmin)
    ax.imshow(spect.T, origin="lower", extent=(kmin,kmax,wmin,wmax), \
        vmin=0.0, vmax=smax, cmap=cm.terrain, interpolation='bilinear', aspect=asp)

    if len(kpoints) > 0:
        knode = np.linspace(0,1,len(kpoints))
        for i in range(1,len(kpoints)-1):
            ax.axvline(x=knode[i],c='k', ls='--')

        ax.set_xticks(knode)        

    if len(klabel) > 0:
        ax.set_xticklabels(klabel)

    ax.set_ylabel(r'$E$ (eV)')

    plt.show()
#----------------------------------------------------------------------
def PlotCDSpectrumPath(ebind,spect,fac=1.0,kps=[],kpoints=[],klabel=[]):
    Ry = 27.211386

    fig, ax = plt.subplots(figsize=(5,5))

    kmin = 0
    kmax = 1
    if len(kps) > 0:
        kmin = np.amin(kps)
        kmax = np.amax(kps)


    smax = fac * np.amax(np.abs(spect))
    wmin, wmax = Ry * np.amin(ebind), Ry * np.amax(ebind)
    asp =  (kmax - kmin) / (wmax - wmin)
    ax.imshow(spect.T, origin="lower", extent=(kmin,kmax,wmin,wmax), \
        vmin=-smax, vmax=smax, cmap=cm.bwr, interpolation='bilinear', aspect=asp)

    if len(kpoints) > 0:
        knode = np.linspace(0,1,len(kpoints))
        for i in range(1,len(kpoints)-1):
            ax.axvline(x=knode[i],c='k', ls='--')

        ax.set_xticks(knode)        

    if len(klabel) > 0:
        ax.set_xticklabels(klabel)

    ax.set_ylabel(r'$E$ (eV)')

    plt.show()

#----------------------------------------------------------------------
def main(argv):
    Ry = 27.211386

    wphot = 20.0 / Ry
    epe, spect = ReadSpectrum(argv[0])

    points = [[-0.5, -0.5, 0.0], [-1./3., -1./3., 0.0], [0.0, 0.0, 0.0], \
        [1./3., 1./3., 0.0], [0.5, 0.5, 0.0]]
    klabel = ["M'", "K'", r"$\Gamma$", "K", "M"]

    ebind = epe - wphot
    PlotSpectrumPath(ebind,spect,fac=1.0,kpoints=points,klabel=klabel)

    if len(argv) > 1:
        epe, spect2 = ReadSpectrum(argv[1])  
        
        ebind = epe - wphot
        PlotSpectrumPath(ebind,spect2,fac=1.0,kpoints=points,klabel=klabel)              

        diff = spect - spect2
        PlotCDSpectrumPath(ebind,diff,fac=1.0,kpoints=points,klabel=klabel)            
#----------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])