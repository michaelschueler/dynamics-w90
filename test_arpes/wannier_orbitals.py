import numpy as np
import h5py
#--------------------------------------------------------------------------------------
class wannorb_slater(object):
	"""docstring for wannorb_slater"""
	#================================================================
	def __init__(self, Ls, Ms, Ns, Zs, Atoms, real_lm=True):
		super(wannorb_slater, self).__init__()
		self.norb = len(Ls)
		if len(Ms) != self.norb:
			print("[Error] len(Ms) != self.norb")
		if len(Ns) != self.norb:
			print("[Error] len(Ns) != self.norb")
		if len(Zs) != self.norb:
			print("[Error] len(Zs) != self.norb")
		if len(Atoms) != self.norb:
			print("[Error] len(Atoms) != self.norb")

		self.l_indx = np.array(Ls)		
		self.m_indx = np.array(Ms)
		self.n_indx = np.array(Ns)
		self.zeff = np.array(Zs)
		self.atom_indx = np.array(Atoms)
		self.real_lm = real_lm
	#================================================================
	def SaveToHDF5(self,fname):
		f = h5py.File(fname, "w")

		f.attrs['wf_type'] = 0
		f.attrs['natoms'] = len(np.unique(self.atom_indx))
		f.attrs['norb'] = self.norb

		f.create_dataset('atom_indx', data=self.atom_indx)
		f.create_dataset('l_indx', data=self.l_indx)
		f.create_dataset('m_indx', data=self.m_indx)
		f.create_dataset('n_indx', data=self.n_indx)
		f.create_dataset('zeff', data=self.zeff)

		if self.real_lm:
			f.attrs['real_lm'] = 1
		else:
			f.attrs['real_lm'] = 0
	
		f.close()
#--------------------------------------------------------------------------------------