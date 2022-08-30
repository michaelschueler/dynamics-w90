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
def run_calc(wphot,mu,Eshift,Epe,Zsc,alpha,phi=0.0):

	kxr = [-1.0, 1.0]
	kyr = [-1.0, 1.0]
	nk1, nk2 = 100, 100

	file_kpts = PrepKpoints(kxr,kyr,nk1,nk2)

	file_ham = "wse2_11bnd_proj_tb.dat"
	file_xyz = "wse2_11bnd_proj.xyz"
	file_orbs = "orbitals_real.h5"
	kpts_reduced = False

	eta = 1.0e-3

	geom = Geometry(alpha,phi=phi)
	pol_circ_p = geom.Get_pol_circ(+1)
	pol_circ_m = geom.Get_pol_circ(-1)

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
		'Zeff': Zsc,
		'polvec_real': list(np.real(pol_circ_p)),
		'polvec_imag': list(np.imag(pol_circ_p))		
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

	file_inp = "inp/wse2_arpes_isocut_pol_circ_p.inp"

	with open(file_inp, 'w') as nml_file:
		f90nml.write(inp, nml_file)

	outpref = "out/wse2_arpes_isocut_pol_circ_p"
	os.system(exe + ' ' + file_inp + ' ' + outpref)

	pes_param['polvec_real'] = list(np.real(pol_circ_m))
	pes_param['polvec_imag'] = list(np.imag(pol_circ_m))	

	inp = {
		'HAMILTONIAN': ham_param,
		'PESPARAMS': pes_param,
		'KPOINTS': kpoints
	}

	file_inp = "inp/wse2_arpes_isocut_pol_circ_m.inp"
	with open(file_inp, 'w') as nml_file:
		f90nml.write(inp, nml_file)

	outpref = "out/wse2_arpes_isocut_pol_circ_m"
	os.system(exe + ' ' + file_inp + ' ' + outpref)	
#---------------------------------------------------------------------
def main():
	Ry = 27.211

	alpha = 0.0 / 180.0 * np.pi
	phi = 0.0 / 180.0 * np.pi

	wphot = 20.0 / Ry
	Eshift = -3.0 / Ry
	mu = -3.0 / Ry

	Emin, Emax = -4.6, -4.2
	Nepe = 4
	Epe = wphot + np.linspace(Emin, Emax, Nepe) / Ry

	Zsc = 4.0
	run_calc(wphot,mu,Eshift,Epe,Zsc,alpha,phi=phi)
#---------------------------------------------------------------------
if __name__ == '__main__':
	main()


