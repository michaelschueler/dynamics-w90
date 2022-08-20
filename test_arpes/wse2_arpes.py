import os
import sys
import numpy as np
import f90nml
from geometry_tools import Geometry
#---------------------------------------------------------------------
def run_calc(wphot,Eshift,Epe,Zsc,alpha,phi=0.0):

	file_ham = "wse2_11bnd_proj_tb.dat"
	file_xyz = "wse2_11bnd_proj.xyz"
	file_orbs = "orbitals_real.h5"
	kpts_reduced = False

	eta = 1.0e-3

	geom = Geometry(alpha,phi=phi)
	pol_p = geom.get_pol_circ(+1)
	pol_m = geom.get_pol_circ(-1)

	ham_param = {
		'file_ham': file_ham,
		'file_xyz': file_xyz
	}

	pes_param = {
		'file_orbs': file_orbs,
		'kpts_reduced': kpts_reduced,
		'gauge': 0,
		'scatt_type': 0,
		'Nepe': len(Epe),
		'wphot': wphot,
		'Eshift': Eshift,
		'Epe_min': np.amin(Epe),
		'Epe_max': np.amax(Epe),
		'lambda_esc': 0.0,
		'eta_smear': eta,
		'Zeff': Zsc,
		'polvec_real': list(np.real(pol_p)),
		'polvec_imag': list(np.imag(pol_p))		
	}

	exe = os.environ['HOME'] + '/Programs/dynamics-w90/exe/arpes.x'

	file_inp = "inp/wse2_arpes_path_polp.inp"
	outpref = "out/wse2_arpes_path_polp"
	os.system(exe + ' ' + file_inp + ' ' + outpref)
#---------------------------------------------------------------------
def main():
	Ry = 27.211

	wphot = 20.0 / Ry
	Eshift = -2.0 / Ry

	Emin, Emax = -6.0, 0.0



