import numpy as np
#--------------------------------------------------------------------------------------
def Rot_z(phi):
    R = np.zeros([3,3])
    R[0,0] = np.cos(phi)
    R[0,1] = np.sin(phi)
    R[1,0] = -np.sin(phi)
    R[1,1] = np.cos(phi)
    R[2,2] = 1.0
    return R
#--------------------------------------------------------------------------------------
class Geometry(object):
    """docstring for rixs_geometry"""
    #================================================================
    def __init__(self, alpha, phi=0.0):
        super(Geometry, self).__init__()
        self.phi = phi
        self.alpha = alpha
        self.Rz = Rot_z(self.phi)
        self.sa = np.sin(self.alpha)
        self.ca = np.cos(self.alpha)
    #================================================================
    def Get_prop_dir(self):
        u = np.array([0.0, -self.sa, -self.ca])
        return self.Rz @ u
    #================================================================
    def Get_pol_p(self):
        e = np.array([0.0, self.ca, -self.sa])
        return self.Rz @ e
    #================================================================
    def Get_pol_s(self):
        e = np.array([1.0, 0.0, 0.0])
        return self.Rz @ e
    #================================================================
    def Get_pol_circ(self,s):
        es = self.Get_pol_s()
        ep = self.Get_pol_p()
        return (es + s * 1j * ep) / np.sqrt(2.0)
    #================================================================
#--------------------------------------------------------------------------------------