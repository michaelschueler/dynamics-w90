import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection
from read_wann_calc_txt import ReadBands, ReadSpin, ReadBerry, ReadOAM
#----------------------------------------------------------------------
def ReadData(fname,nk,nbnd):
    data = np.loadtxt(fname)

    epsk = np.zeros([nk,nbnd])
    oam = np.zeros([nk,nbnd])

    for ik in range(nk):
        epsk[ik,:] = data[ik*nbnd:(ik+1)*nbnd,2]
        oam[ik,:] = data[ik*nbnd:(ik+1)*nbnd,3]

    return epsk, oam
#----------------------------------------------------------------------
def PlotBandsOAM(epsk,oam,klabel,knodes,Eshift=0.0,title="",fout=""):
    nk = epsk.shape[0]
    nbnd = epsk.shape[1]

    xk = np.linspace(0, 1, nk)

    fig, ax = plt.subplots()

    for i in range(nbnd):
        points = np.array([xk, epsk[:,i] + Eshift]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        norm = plt.Normalize(-4,4)
        lc = LineCollection(segments, cmap='seismic', norm=norm)
        # Set the values used for colormapping
        lc.set_array(oam[:,i])
        lc.set_linewidth(2)
        line = ax.add_collection(lc)

    ax.set_xlim(0.0,1.0)
    ax.set_ylim(-5, 3)
    ax.set_ylabel(r'$E$ (eV)')

    ax.set_xticks(knodes)
    ax.set_xticklabels(klabel)

    for i in range(1,len(klabel)-1):
        ax.axvline(x=knodes[i],c='k', ls='--')

    if len(title) > 0:
        plt.title(title,fontsize=16)

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#----------------------------------------------------------------------
def PlotOAM(Lz_ref,Lz_calc,klabel,knodes,fout=""):
    nk = Lz_ref.shape[0]
    xk = np.linspace(0,1,nk)

    fig, ax = plt.subplots()

    ax.plot(xk, Lz_ref, c='tab:red', label="ref.")
    ax.plot(xk, Lz_calc, c='tab:blue', label="wannier_calc")

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$L_z \,(\hbar)$')

    for i in range(1,len(klabel)-1):
        ax.axvline(x=knodes[i],c='k', ls='--')

    ax.set_xticks(knodes)
    ax.set_xticklabels(klabel)

    ax.legend(loc='best',frameon=False,fontsize=14)

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#----------------------------------------------------------------------
def main(argv):
    Ry = 27.211386

    fname = "data/WS2_soc/orbital.dat"

    nk, nbnd = 121, 80
    epsk, oam = ReadData(fname,nk,nbnd)
    xk = np.linspace(0, 1, nk)

    Evbm = -1.4875
    klabel = [r"$\Gamma$", "K", "M", r"$\bar{\mathrm{K}}$", r"$\Gamma$"]
    knodes = [xk[0], xk[40], xk[60], xk[80], xk[-1]]
    PlotBandsOAM(epsk,oam,klabel,knodes,Eshift=-Evbm,title="ref.")

    ibnd = 25
    Lz_ref = oam[:,ibnd]

    fname = "out/ws2_soc"
    epsk = ReadBands(fname)
    oam = ReadOAM(fname)

    Evbm = -1.4875
    klabel = [r"$\Gamma$", "K", "M", r"$\bar{\mathrm{K}}$", r"$\Gamma$"]
    PlotBandsOAM(Ry*epsk,oam[:,2,:],klabel,knodes,Eshift=-Evbm,title="wannier_calc")

    ibnd = 13
    Lz_calc = oam[:,2,ibnd]

    PlotOAM(Lz_ref,Lz_calc,klabel,knodes)
#----------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])