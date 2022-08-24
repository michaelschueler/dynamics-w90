import os
import sys
import numpy as np
from wannier_orbitals import wannorb_slater
#--------------------------------------------------------------------------------------
def main():
    Ry = 27.211386

    # set up Wannier orbitals
    #    dz2  dxz  dyz  dx2y2  dxy  pz  px  py  pz  px  py                                                                                     
    Ls = [ 2,   2,   2,     2,   2,  1,  1,  1,  1,  1,  1]
    Ms = [ 0,   1,  -1,     2,  -2,  0,  1, -1,  0,  1, -1]
    Ns = [ 3,   3,   3,     3,   3,  2,  2 , 2,  2,  2,  2]
    Zorb = np.ones(len(Ls))
    Zscatt = np.concatenate((0.2*np.ones(5), 0.2*np.ones(6)))
    Atoms = [74, 74, 74, 74, 74, 34, 34, 34, 34, 34, 34]

    orb = wannorb_slater(Ls, Ms, Ns, Zorb, Atoms, Zscatt=Zscatt, real_lm=True)

    file_wannorb = "orbitals_real.h5"
    orb.SaveToHDF5(file_wannorb)

    weight = np.concatenate((np.ones(5), np.zeros(6)))
    orb = wannorb_slater(Ls, Ms, Ns, Zorb, Atoms, weight=weight, Zscatt=Zscatt, real_lm=True)
    file_wannorb = "orbitals_donly_real.h5"
    orb.SaveToHDF5(file_wannorb)
#--------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
    
    



