import sys
import numpy as np
import h5py
#----------------------------------------------------------------------
def ReadParams(fname):
    f = h5py.File(fname, "r")
    nk = f.attrs['nk']
    nwan = f.attrs['nwan']
    nbnd = f.attrs['nbnd']
    f.close()
    return nk, nwan, nbnd
#----------------------------------------------------------------------
def ReadBands(fname):
    f = h5py.File(fname, "r")
    epsk = np.array(f['epsk'])
    f.close()
    return epsk
#----------------------------------------------------------------------
def ReadKpts(fname):
    f = h5py.File(fname, "r")
    kpts = np.array(f['kpts'])
    f.close()
    return kpts
#----------------------------------------------------------------------
def ReadOrbWeight(fname):
    f = h5py.File(fname, "r")
    orb_weight = np.array(f['orbweight'])
    f.close()
    return orb_weight
#----------------------------------------------------------------------
def ReadSpin(fname):
    f = h5py.File(fname, "r")
    spin = np.array(f['spin'])
    f.close()
    return spin
#----------------------------------------------------------------------
def ReadBerry(fname):
    f = h5py.File(fname, "r")
    berry = np.array(f['berry'])
    f.close()
    return berry
#----------------------------------------------------------------------
def ReadOAM(fname):
    f = h5py.File(fname, "r")
    oam = np.array(f['oam'])
    f.close()
    return oam
#----------------------------------------------------------------------
def ReadMetric(fname):
    f = h5py.File(fname, "r")
    oam = np.array(f['metric'])
    f.close()
    return oam
#----------------------------------------------------------------------
def ReadEvecs(fname):
    f = h5py.File(fname, "r")
    evecs = np.array(f['evecs-real']) + 1j * np.array(f['evecs-imag'])
    f.close()
    return evecs
#----------------------------------------------------------------------
