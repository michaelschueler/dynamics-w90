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
def PlotSpectrum(spect,kxr,kyr,nk1,nk2,icut=0,fac=1.0):
    Ry = 27.211386

    fig, ax = plt.subplots(figsize=(5,5))

    smax = fac * np.amax(spect)

    data = np.reshape(spect[:,icut], [nk1,nk2])
    ax.imshow(data.T, origin="lower", extent=(kxr[0],kxr[1],kyr[0],kyr[1]), \
        vmin=0.0, vmax=smax, cmap=cm.terrain, interpolation='bilinear')

    plt.show()
#----------------------------------------------------------------------
def PlotCDSpectrum(spect,kxr,kyr,nk1,nk2,icut=0,fac=1.0):
    Ry = 27.211386

    fig, ax = plt.subplots(figsize=(5,5))

    smax = fac * np.amax(np.abs(spect))

    data = np.reshape(spect[:,icut], [nk1,nk2])
    ax.imshow(data.T, origin="lower", extent=(kxr[0],kxr[1],kyr[0],kyr[1]), \
        vmin=-smax, vmax=smax, cmap=cm.bwr, interpolation='bilinear')

    plt.show()

#----------------------------------------------------------------------
def main(argv):
    Ry = 27.211386

    kxr = [-1.0, 1.0]
    kyr = [-1.0, 1.0]
    nk1, nk2 = 100, 100

    wphot = 20.0 / Ry
    epe, spect = ReadSpectrum(argv[0])

    for icut in range(spect.shape[1]):
        PlotSpectrum(spect,kxr,kyr,nk1,nk2,fac=1.0,icut=icut)

    if len(argv) > 1:
        epe, spect2 = ReadSpectrum(argv[1])  
        
        for icut in range(spect.shape[1]):
            PlotSpectrum(spect2,kxr,kyr,nk1,nk2,fac=1.0,icut=icut)

        diff = spect - spect2
        for icut in range(spect.shape[1]):
            PlotCDSpectrum(diff,kxr,kyr,nk1,nk2,fac=0.2,icut=icut)            
#----------------------------------------------------------------------
if __name__ == '__main__':
    main(sys.argv[1:])