import sys
import numpy as np
import matplotlib.pyplot as plt
from read_wann_calc import ReadOAM, ReadEvecs
#----------------------------------------------------------------------
def GetMRot():
    isq2 = 1.0/np.sqrt(2.0)
    rot = np.zeros([11, 11], dtype=np.complex_)
    rot[0,2] = 1.0
    rot[1,1] = isq2
    rot[1,3] = -isq2
    rot[2,1] = 1j * isq2
    rot[2,3] = 1j * isq2
    rot[3,0] = isq2
    rot[3,4] = isq2
    rot[4,0] = 1j * isq2
    rot[4,4] = -1j * isq2

    for i in range(5,11):
        rot[i,i] = 1.0

    return rot
#----------------------------------------------------------------------
def GetOAM_local(evecs):

    mrot = GetMRot()
    evecs_mb = np.einsum('rm,kar->kma', mrot, evecs)

    ms = np.arange(-2, 3, step=1)
    Lz = np.einsum('m, kma -> ka', ms, np.abs(evecs_mb[:,0:5,:])**2)

    return Lz
#----------------------------------------------------------------------
def PlotOAM(Lz_berry,Lz_loc,ibnd,klabel,fout=""):
    nk = Lz_berry.shape[0]
    xk = np.linspace(0,1,nk)

    fig, ax = plt.subplots()

    ax.plot(xk, Lz_berry[:,ibnd], c='tab:red', label="Berry")
    ax.plot(xk, Lz_loc[:,ibnd], c='tab:blue', label="local (W)")

    ax.set_xlim(0.0,1.0)
    ax.set_ylabel(r'$L_z \,(\hbar)$')

    knode = np.linspace(0,1,len(klabel))
    for i in range(1,len(klabel)-1):
        ax.axvline(x=knode[i],c='k', ls='--')

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    ax.legend(loc='best',frameon=False,fontsize=14)

    if len(fout) > 0:
        plt.savefig(fout, bbox_inches="tight")
    else:
        plt.show()
#----------------------------------------------------------------------
def main(argv):
    if len(argv) > 0:
        fname = argv[0]
    else:
        print("No input. Exiting ...")
        exit()

    fout = ""
    if len(argv) > 1:
        fout = argv[1]

    klabel = ["M'", "K'", r"$\Gamma$", "K", "M"]

    evecs = ReadEvecs(fname)
    Lz_loc = GetOAM_local(evecs)

    oam = ReadOAM(fname)
    Lz_berry = oam[:,2,:]

    ibnd = 6
    PlotOAM(Lz_berry,Lz_loc,ibnd,klabel,fout=fout)
#----------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])