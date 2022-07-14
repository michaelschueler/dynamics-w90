import numpy as np
#--------------------------------------------------------------------------------------
def GenLmatrix_mbasis(L0):
	nl = 2*L0+1

	Lmat = np.zeros([nl,nl,3],dtype=np.complex_)

	ms = np.arange(-L0, L0+1, step=1)
	Lp = np.zeros([nl, nl])
	Lm = np.zeros([nl, nl])
	for i,m in enumerate(ms[0:nl-1]):
		Lp[i,i+1] = np.sqrt((L0 - m)*(L0 + m + 1))
	for i,m in enumerate(ms[1:]):
		Lm[i+1,i] = np.sqrt((L0 + m) * (L0 - m + 1))

	Lmat[:,:,0] = 0.5*(Lp + Lm)
	Lmat[:,:,1] = -0.5j*(Lp - Lm)
	Lmat[:,:,2] = np.diag(ms)

	return Lmat
#--------------------------------------------------------------------------------------
def Transform_w90(L0,Lmat):
	nl = 2*L0+1
	# Lmat_w90 = np.zeros([nl,nl,3], dtype=np.complex_)

	Arot = np.zeros([nl,nl], dtype=np.complex_)

	sq2 = np.sqrt(2.0)

	if L0 == 0:
		Arot[0,0] = 1.0
	elif L0 == 1:
		Arot[0,0+L0] = 1.0
		Arot[1,-1+L0] = 1.0/sq2
		Arot[1,1+L0] = -1.0/sq2
		Arot[2,-1+L0] = 1j/sq2
		Arot[2,1+L0] = 1j/sq2
	elif L0 == 2:
		Arot[0,0+L0] = 1.0
		Arot[1,-1+L0] = 1.0/sq2
		Arot[1,1+L0] = -1.0/sq2
		Arot[2,-1+L0] = 1j/sq2
		Arot[2,1+L0] = 1j/sq2
		Arot[3,-2+L0] = 1.0/sq2
		Arot[3,+2+L0] = 1.0/sq2
		Arot[4,-2+L0] = 1j/sq2
		Arot[4,+2+L0] = -1j/sq2
	elif L0 == 3:
		Arot[0,0+L0] = 1.0
		Arot[1,-1+L0] = 1.0/sq2
		Arot[1,+1+L0] = -1.0/sq2
		Arot[2,-1+L0] = 1j/sq2
		Arot[2,+1+L0] = 1j/sq2
		Arot[3,-2+L0] = 1.0/sq2
		Arot[3,+2+L0] = 1.0/sq2
		Arot[4,-2+L0] = 1j/sq2
		Arot[4,+2+L0] = -1j/sq2
		Arot[5,-3+L0] = 1.0/sq2
		Arot[5,+3+L0] = -1.0/sq2
		Arot[6,-3+L0] = 1j/sq2
		Arot[6,+3+L0] = 1j/sq2
	else:
		print("[Transform_w90] higher angular momenta not implemented yet!")
		exit()

	Lmat_w90 = np.einsum('im,mnr,jn->ijr',Arot,Lmat,np.conj(Arot))

	# print(Lmat_w90[:,:,0])
	# print(Lmat_w90[:,:,1])
	# print(Lmat_w90[:,:,2])
	# exit()

	return Lmat_w90
#--------------------------------------------------------------------------------------
def GenLmatrix_w90basis(L0):

	Lmat_mbasis = GenLmatrix_mbasis(L0)
	Lmat_w90 = Transform_w90(L0,Lmat_mbasis)
	return Lmat_w90
#--------------------------------------------------------------------------------------
def WriteSOCData(fname,Ls,Lmats):
	ngroups = len(Ls)
	if len(Lmats) != ngroups:
		print("[WriteSOCData] len(Lmats) != ngroups")
		exit()

	ndim = 2*np.array(Ls) + 1
	norb = np.sum(ndim)
	hsoc = np.zeros([norb,norb,3], dtype=np.complex_)

	orb_indx = np.zeros(ngroups,dtype=np.int_)
	orb_indx[1:] = np.cumsum(ndim)[0:-1]

	for i in range(ngroups):
		hsoc[orb_indx[i]:orb_indx[i]+ndim[i],orb_indx[i]:orb_indx[i]+ndim[i],:] = Lmats[i]

	f = open(fname, "w")
	f.write("{}\n".format(ngroups))
	for i in range(ngroups):
		f.write("{}  ".format(ndim[i]))
	f.write("\n")
	for i in range(ngroups):
		f.write("{}  ".format(Ls[i]))
	f.write("\n")

	for i in range(norb):
		for j in range(norb):
			f.write("{:14.7f}  {:14.7f}  {:14.7f}  {:14.7f}  {:14.7f}  {:14.7f}\n".\
				format(np.real(hsoc[i,j,0]), np.real(hsoc[i,j,1]), np.real(hsoc[i,j,2]), \
				 np.imag(hsoc[i,j,0]), np.imag(hsoc[i,j,1]), np.imag(hsoc[i,j,2])))
	f.close()
#--------------------------------------------------------------------------------------
def WriteSOCLam(fname,lam):
	f = open(fname, "w")
	f.write(str(len(lam)) + "\n")
	for i in range(len(lam)):
		f.write("{:14.7f} ".format(lam[i]))
	f.close()
#--------------------------------------------------------------------------------------
if __name__ == '__main__':
	Ls = [1, 2]
	Lmat1 = GenLmatrix_mbasis(Ls[0])
	Lmat2 = GenLmatrix_mbasis(Ls[1])

	Lmats = [Lmat1, Lmat2]

	WriteSOCData("", Ls, Lmats)

