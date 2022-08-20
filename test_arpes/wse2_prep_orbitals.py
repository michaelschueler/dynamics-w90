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
    Zs = np.ones(len(Ls))
    Atoms = [74, 74, 74, 74, 74, 34, 34, 34, 34, 34, 34]
    orb = wannorb_slater(Ls, Ms, Ns, Zs, Atoms, real_lm=True)

    file_wannorb = "orbitals_real.h5"
    orb.SaveToHDF5(file_wannorb)

    orb = wannorb_slater(Ls, Ms, Ns, Zs, Atoms, real_lm=False)
    
    file_wannorb = "orbitals_cplx.h5"
    orb.SaveToHDF5(file_wannorb)
#--------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
    
    



