import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
#---------------------------------------------------------------------- 
def Plot_bandstructure(epsk,klabel,fout=""):
    Ry = 27.211386

    nk = epsk.shape[0]
    nbnd = epsk.shape[1]
    xk = np.linspace(0.0, 1.0, nk)

    fig, ax = plt.subplots()

    for ibnd in range(nbnd):
        ax.plot(xk, Ry * epsk[:,ibnd], c='blue')

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(klabel))
    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#---------------------------------------------------------------------- 
def Plot_bandstructure_orbweight(epsk,orb_weight,klabel,iorb=0,Emin=-100.0,Emax=100.0,fout=""):
    Ry = 27.211386

    nk = epsk.shape[0]
    nbnd = epsk.shape[1]
    xk = np.linspace(0.0, 1.0, nk)

    wkorb = orb_weight[:,:,iorb]

    fig, ax = plt.subplots()

    for i in range(epsk.shape[1]):
        points = np.array([xk, Ry*epsk[:,i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(0, 1)
        lc = LineCollection(segments, cmap=cm.afmhot_r, norm=norm)
        # Set the values used for colormapping
        lc.set_array(wkorb[:,i])
        lc.set_linewidth(4)
        line = ax.add_collection(lc)

    plt.colorbar(line)

    Erange = np.amax(epsk) - np.amin(epsk)
    Emin_ = Ry*(np.amin(epsk) - 0.02 * Erange)
    Emax_ = Ry*(np.amax(epsk) + 0.02 * Erange)
    if Emin > -90.0:
        Emin_ = Emin
    if Emax < 90.0:
        Emax_ = Emax
    ax.set_ylim(Emin_,Emax_)

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(klabel))
    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#---------------------------------------------------------------------- 
def Plot_bandstructure_spin(epsk,spin,klabel,spin_dir=2,Emin=-100.0,Emax=100.0,fout=""):
    Ry = 27.211386

    nk = epsk.shape[0]
    nbnd = epsk.shape[1]
    xk = np.linspace(0.0, 1.0, nk)

    Sp = spin[:,:,spin_dir]

    fig, ax = plt.subplots()

    for i in range(epsk.shape[1]):
        points = np.array([xk, Ry*epsk[:,i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(-1, 1)
        lc = LineCollection(segments, cmap='bwr', norm=norm)
        # Set the values used for colormapping
        lc.set_array(Sp[:,i])
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

    plt.colorbar(line)

    Erange = np.amax(epsk) - np.amin(epsk)
    Emin_ = Ry*(np.amin(epsk) - 0.02 * Erange)
    Emax_ = Ry*(np.amax(epsk) + 0.02 * Erange)
    if Emin > -90.0:
        Emin_ = Emin
    if Emax < 90.0:
        Emax_ = Emax
    ax.set_ylim(Emin_,Emax_)

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(klabel))

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#---------------------------------------------------------------------- 
def Plot_bandstructure_berry(epsk,berry,klabel,idir=2,Bmax=1.0,Emin=-100.0,Emax=100.0,fout=""):
    Ry = 27.211386

    nk = epsk.shape[0]
    nbnd = epsk.shape[1]
    xk = np.linspace(0.0, 1.0, nk)

    wb = berry[:,idir,:]

    fig, ax = plt.subplots()

    for i in range(epsk.shape[1]):
        points = np.array([xk, Ry*epsk[:,i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(-Bmax, Bmax)
        lc = LineCollection(segments, cmap='bwr', norm=norm)
        # Set the values used for colormapping
        lc.set_array(wb[:,i])
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

    plt.colorbar(line)

    Erange = np.amax(epsk) - np.amin(epsk)
    Emin_ = Ry*(np.amin(epsk) - 0.02 * Erange)
    Emax_ = Ry*(np.amax(epsk) + 0.02 * Erange)
    if Emin > -90.0:
        Emin_ = Emin
    if Emax < 90.0:
        Emax_ = Emax
    ax.set_ylim(Emin_,Emax_)

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(klabel))

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#---------------------------------------------------------------------- 
def Plot_bandstructure_oam(epsk,oam,klabel,idir=2,Lmax=1.0,Emin=-100.0,Emax=100.0,fout=""):
    Ry = 27.211386

    nk = epsk.shape[0]
    nbnd = epsk.shape[1]
    xk = np.linspace(0.0, 1.0, nk)

    Lb = oam[:,idir,:]

    fig, ax = plt.subplots()

    for i in range(epsk.shape[1]):
        points = np.array([xk, Ry*epsk[:,i]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(-Lmax, Lmax)
        lc = LineCollection(segments, cmap='seismic', norm=norm)
        # Set the values used for colormapping
        lc.set_array(Lb[:,i])
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

    plt.colorbar(line)

    Erange = np.amax(epsk) - np.amin(epsk)
    Emin_ = Ry*(np.amin(epsk) - 0.02 * Erange)
    Emax_ = Ry*(np.amax(epsk) + 0.02 * Erange)
    if Emin > -90.0:
        Emin_ = Emin
    if Emax < 90.0:
        Emax_ = Emax
    ax.set_ylim(Emin_,Emax_)

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(klabel))

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#---------------------------------------------------------------------- 