import sys
import numpy as np
#----------------------------------------------------------------------
def ReadBands(pref):
    epsk = np.loadtxt(pref + "_epsk.txt")
    return epsk
#----------------------------------------------------------------------
def ReadKpts(pref):
    kpts = np.loadtxt(pref + "_kpts.txt")
    return kpts
#----------------------------------------------------------------------
def ReadOrbWeight(pref,norb):
    w = np.loadtxt(pref + "_orbweight_1.txt")
    nk = w.shape[0]
    nbnd = w.shape[1]
    worb = np.zeros([nk,nbnd,norb])
    for i in range(norb):
        fname = pref + "_orbweight_{}.txt".format(i+1)
        w = np.loadtxt(fname)
        worb[:,:,i] = w
    return worb
#----------------------------------------------------------------------
def ReadSpin(pref):
    spin_x = np.loadtxt(pref + "_spin_x.txt")
    spin_y = np.loadtxt(pref + "_spin_y.txt")
    spin_z = np.loadtxt(pref + "_spin_z.txt")

    nk = spin_x.shape[0]
    nbnd = spin_x.shape[1]
    spin = np.zeros([nk,nbnd,3])
    spin[:,:,0] = spin_x
    spin[:,:,1] = spin_y
    spin[:,:,2] = spin_z

    return spin
#----------------------------------------------------------------------
def ReadBerry(pref):
    berry_x = np.loadtxt(pref + "_berry_x.txt")
    berry_y = np.loadtxt(pref + "_berry_y.txt")
    berry_z = np.loadtxt(pref + "_berry_z.txt")

    nk = berry_x.shape[0]
    nbnd = berry_x.shape[1]
    berry = np.zeros([nk,3,nbnd])
    berry[:,0,:] = berry_x
    berry[:,1,:] = berry_y
    berry[:,2,:] = berry_z

    return berry
#----------------------------------------------------------------------
def ReadOAM(pref):
    oam_x = np.loadtxt(pref + "_oam_x.txt")
    oam_y = np.loadtxt(pref + "_oam_y.txt")
    oam_z = np.loadtxt(pref + "_oam_z.txt")

    nk = oam_x.shape[0]
    nbnd = oam_x.shape[1]
    oam = np.zeros([nk,3,nbnd])
    oam[:,0,:] = oam_x
    oam[:,1,:] = oam_y
    oam[:,2,:] = oam_z

    return oam
#----------------------------------------------------------------------
