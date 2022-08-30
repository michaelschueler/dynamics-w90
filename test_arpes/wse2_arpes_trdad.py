import os
import sys
import numpy as np
import f90nml
from geometry_tools import Geometry
#---------------------------------------------------------------------
def PrepKpoints(kxr,kyr,nk1,nk2):
	file_kpts = "inp/kpatch_x{}_{}_y{}_{}_nk{}x{}.dat"\
		.format(kxr[0],kxr[1],kyr[0],kyr[1],nk1,nk2)

	kxs = np.linspace(kxr[0], kxr[1], nk1)
	kys = np.linspace(kyr[0], kyr[1], nk2)

	X, Y = np.meshgrid(kxs, kys)
	xx = np.reshape(X, [nk1*nk2])
	yy = np.reshape(Y, [nk1*nk2])	

	np.savetxt(file_kpts, np.c_[xx, yy])
	
	return file_kpts
#---------------------------------------------------------------------
def run_calc(wphot,mu,Eshift,Epe,alpha,phi=0.0):

	kxr = [-1.0, 1.0]
	kyr = [-1.0, 1.0]
	nk1, nk2 = 100, 100

	file_kpts = PrepKpoints(kxr,kyr,nk1,nk2)

	file_ham = "wse2_11bnd_proj_tb.dat"
	file_xyz = "wse2_11bnd_proj.xyz"
	file_orbs = "orbitals_real.h5"
	kpts_reduced = False

	eta = 1.0e-3

	geom = Geometry(np.pi*alpha/180,phi=np.pi*phi/180)
	pol_p = geom.Get_pol_p()

	ham_param = {
		'file_ham': file_ham,
		'file_xyz': file_xyz,
		'MuChem': mu
	}

	pes_param = {
		'file_orbs': file_orbs,
		'kpts_reduced': kpts_reduced,
		'gauge': 1,
		'scatt_type': 1,
		'Nepe': len(Epe),
		'wphot': wphot,
		'Eshift': Eshift,
		'Epe_min': np.amin(Epe),
		'Epe_max': np.amax(Epe),
		'lambda_esc': 0.0,
		'eta_smear': eta,
		'polvec_real': list(np.real(pol_p)),
		'polvec_imag': list(np.imag(pol_p))		
	}

	kpoints = {
		'kpoints_type': 'list',
		'file_kpts': file_kpts
	}

	inp = {
		'HAMILTONIAN': ham_param,
		'PESPARAMS': pes_param,
		'KPOINTS': kpoints
	}

	exe = os.environ['HOME'] + '/Programs/dynamics-w90/exe/arpes.x'

	file_inp = "inp/wse2_arpes_isocut_alpha{}_phi{}_pol_p.inp".format(alpha,phi)
	with open(file_inp, 'w') as nml_file:
		f90nml.write(inp, nml_file)

	outpref = "out/wse2_arpes_isocut_alpha{}_phi{}_pol_p".format(alpha,phi)
	os.system(exe + ' ' + file_inp + ' ' + outpref)
#---------------------------------------------------------------------
def main():
	Ry = 27.211

	alpha = 65
	phi = 30

	wphot = 20.0 / Ry
	Eshift = -3.0 / Ry
	mu = -3.0 / Ry

	Emin, Emax = -4.6, -4.2
	Nepe = 4
	Epe = wphot + np.linspace(Emin, Emax, Nepe) / Ry

	Zsc = 1.5
	run_calc(wphot,mu,Eshift,Epe,alpha,phi=phi)

	alpha = -65
	phi = -30
	run_calc(wphot,mu,Eshift,Epe,alpha,phi=phi)
#---------------------------------------------------------------------
if __name__ == '__main__':
	main()


