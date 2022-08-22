import os
home = os.environ['HOME']
import sys
sys.path.append(home + '/Programs/pythonw90/')
import numpy as np
import matplotlib.pyplot as plt
from kgrid import GenKpath
from python_w90 import wann90_tb
#----------------------------------------------------------------------
def PlotWanBands(file_ham,points,subplot,Eshift=0.0,col='blue',nseg=80):
    Ry = 27.211386

    wann = wann90_tb(file_ham)
    kpts = GenKpath(points,nseg)
    epsk = wann.get_energy(kpts)
    
    Nk = len(epsk)
    xk = np.linspace(0.0,1.0,Nk)

    for ibnd in range(epsk.shape[1]):
        subplot.plot(xk,Ry*epsk[:,ibnd] + Eshift, c=col)
#----------------------------------------------------------------------
def main():
    Ry = 27.211386

    file_wan = "wse2_11bnd_proj_tb.dat"

    fig, ax = plt.subplots(figsize=(5,5))

    points = [[-0.5, -0.5, 0.0], [-1./3., -1./3., 0.0], [0.0, 0.0, 0.0], \
        [1./3., 1./3., 0.0], [0.5, 0.5, 0.0]]
    klabel = ["M'", "K'", r"$\Gamma$", "K", "M"]

    Eshift = -3.0

    PlotWanBands(file_wan,points,ax,Eshift=Eshift,col='blue')

    ax.set_xlim(0.0,1.0)
    #ax.set_ylim(-2.4,-2.0)
    ax.set_ylabel(r'$E$ (eV)')

    knode = np.linspace(0,1,len(points))
    for i in range(1,len(points)-1):
            ax.axvline(x=knode[i],c='k', ls='--')

    ax.set_xticks(knode)
    ax.set_xticklabels(klabel)

    plt.show()
#----------------------------------------------------------------------
if __name__ == '__main__':
    main()