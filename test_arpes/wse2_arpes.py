import os
import sys
import numpy as np
import f90nml
from geometry_tools import Geometry
#---------------------------------------------------------------------
def run_calc(wphot,mu,Eshift,Epe,Zsc,alpha,phi=0.0):

	file_ham = "wse2_11bnd_proj_tb.dat"
	file_xyz = "wse2_11bnd_proj.xyz"
	file_orbs = "orbitals_real.h5"
	file_kpts = "path_MKGKM.inp"
	kpts_reduced = True

	eta = 4.0e-3

	geom = Geometry(alpha,phi=phi)
	pol_p = geom.Get_pol_circ(+1)
	pol_m = geom.Get_pol_circ(-1)

	ham_param = {
		'file_ham': file_ham,
		'file_xyz': file_xyz,
		'MuChem': mu
	}

	pes_param = {
		'file_orbs': file_orbs,
		'kpts_reduced': kpts_reduced,
		'lambda_orbital_term': False,
		'gauge': 1,
		'scatt_type': 0,
		'Nepe': len(Epe),
		'wphot': wphot,
		'Eshift': Eshift,
		'Epe_min': np.amin(Epe),
		'Epe_max': np.amax(Epe),
		'lambda_esc': 0.0,
		'eta_smear': eta,
		'radint_numpoints_k': 10,
		'radint_numpoints_r': 128,
		'polvec_real': list(np.real(pol_p)),
		'polvec_imag': list(np.imag(pol_p))		
	}

	kpoints = {
		'kpoints_type': 'path',
		'file_kpts': file_kpts
	}

	inp = {
		'HAMILTONIAN': ham_param,
		'PESPARAMS': pes_param,
		'KPOINTS': kpoints
	}

	exe = os.environ['HOME'] + '/Programs/dynamics-w90/exe/arpes.x'

	file_inp = "inp/wse2_arpes_path_polp.inp"

	with open(file_inp, 'w') as nml_file:
		f90nml.write(inp, nml_file)

	outpref = "out/wse2_arpes_path_polp"
	os.system(exe + ' ' + file_inp + ' ' + outpref)

	exit()

	pes_param['polvec_real'] = list(np.real(pol_m))
	pes_param['polvec_imag'] = list(np.imag(pol_m))	

	inp = {
		'HAMILTONIAN': ham_param,
		'PESPARAMS': pes_param,
		'KPOINTS': kpoints
	}

	file_inp = "inp/wse2_arpes_path_polm.inp"
	with open(file_inp, 'w') as nml_file:
		f90nml.write(inp, nml_file)

	outpref = "out/wse2_arpes_path_polm"
	os.system(exe + ' ' + file_inp + ' ' + outpref)	
#---------------------------------------------------------------------
def main():
	Ry = 27.211

	alpha = 0.0 / 180.0 * np.pi

	wphot = 20.0 / Ry
	Eshift = -3.0 / Ry
	mu = -3.0 / Ry

	Emin, Emax = -5.0, -3.0
	Nepe = 20
	Epe = wphot + np.linspace(Emin, Emax, Nepe) / Ry

	Zsc = 0.0
	run_calc(wphot,mu,Eshift,Epe,Zsc,alpha,phi=0.0)
#---------------------------------------------------------------------
if __name__ == '__main__':
	main()


