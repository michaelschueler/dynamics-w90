import sys
import os
import numpy as np
import subprocess
import f90nml
#--------------------------------------------------------------------------------------
class WannierARPES():
    #========================================
    def __init__(self, PathExe, PathInp='./inp/', PathOut='./out/', PathLog='./log/', mpicmd=""):
        self.exe = os.path.join(PathExe, "wann_arpes_mpi.x")
        self.mpicmd = mpicmd
        self.PathInp = PathInp
        self.PathOut = PathOut
        self.PathLog = PathLog
    #========================================
    def SetHamiltonian(self, file_ham:str, MuChem:float=0.0, file_ovlp:str="", 
                       slab_mode:bool=False, slab_nlayers:int=1):
        
        self.ham_param = {
            'file_ham': file_ham,
            'MuChem': MuChem,
            'file_ovlp': file_ovlp,
            'slab_mode': slab_mode
        }
        
        self.slab_param = {}
        if slab_mode:
            self.slab_param['slab_nlayers'] = slab_nlayers
    #========================================
    def SetPESParams(self, params, photon):
        
        file_orbs = params.get('file_orbs', "")
        if len(file_orbs) == 0:
            print("[Error] file_orbs is not set.")
            sys.exit(1)
        
        file_scatt = params.get('file_scatt', "")
        kpts_reduced = params.get('kpts_reduced', False)
        lambda_orbital_term = params.get('lambda_orbital_term', False)
        gauge = params.get('gauge', 'length').strip().lower()
        if gauge not in ['length', 'velocity']:
            print("[Error] gauge must be 'length' or 'velocity'.")
            sys.exit(1)

        scatt_flag = params.get('scatt_type', 'pw').strip().lower()
        if scatt_flag == 'pw':
            scatt_type = 0
        elif scatt_type == 'coulomb':
            scatt_type = 1
        elif scatt_type == 'input':
            scatt_type = 2
        else:
            print("[Error] scatt_type must be 'pw', 'coulomb', or 'input'.")
            sys.exit(1)

        Eshift = params.get('Eshift', 0.0)
        radint_numpoints_k = params.get('radint_numpoints_k', 40)
        radint_numpoints_r = params.get('radint_numpoints_r', 128)

        Nepe = photon.get('Nepe', 100)
        wphot = photon.get('wphot', 1.0)
        Epe_min = photon.get('Epe_min', -1.0)
        Epe_max = photon.get('Epe_max', 1.0)
        lambda_esc = photon.get('lambda_esc', 0.0)
        eta_smear = photon.get('eta_smear', 0.01)
        polvec_real = photon.get('polvec_real', [1.0, 0.0, 0.0])
        polvec_imag = photon.get('polvec_imag', [0.0, 0.0, 0.0])

        self.pes_param = {
            'file_orbs': file_orbs,
            'file_scatt': file_scatt,
            'kpts_reduced': kpts_reduced,
            'lambda_orbital_term': lambda_orbital_term,
            'gauge': gauge,
            'scatt_type': scatt_type,
            'Eshift': Eshift,
            'radint_numpoints_k': radint_numpoints_k,
            'radint_numpoints_r': radint_numpoints_r,
            'Nepe': Nepe,
            'wphot': wphot,
            'Epe_min': Epe_min,
            'Epe_max': Epe_max,
            'lambda_esc': lambda_esc,
            'eta_smear': eta_smear,
            'polvec_real': polvec_real,
            'polvec_imag': polvec_imag
        }
    #========================================
    def SetKPTS(self, kpoints_type:str, file_kpts:str="", nk1:int=1, nk2:int=1):
        self.kpoints = {
            'kpoints_type': kpoints_type,
            'file_kpts': file_kpts,
            'nk1': nk1,
            'nk2': nk2
        }
    #========================================
    def __WriteInput(self, file_inp:str):

        inp = {
            'HAMILTONIAN': self.ham_param,
            'SLAB': self.slab_param,
            'PESPARAMS': self.pes_param,
            'KPOINTS': self.kpoints
        }

        with open(file_inp, 'w') as nml_file:
            f90nml.write(inp, nml_file)
    #========================================
    def Run(self, prefix:str, debug_mode:bool=False):
        
        file_inp = os.path.join(self.PathInp, prefix + '.inp')
        file_out = os.path.join(self.PathOut, prefix)
        file_log = os.path.join(self.PathLog, prefix + '.log')

        self.__WriteInput(file_inp)
        
        if len(self.mpicmd):
            cmd = [self.mpicmd, self.exe, file_inp, file_out]
        else:
            cmd = [self.exe, file_inp, file_out]

        if debug_mode:
            subprocess.run(cmd)
        else:
            with open(file_log, 'w') as log_file:
                subprocess.run(cmd, stdout=log_file, stderr=subprocess.STERR)
        
        if debug_mode:
            subprocess.run(cmd, check=True)
        else:
            with open(file_log, 'w') as log_file:
                subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True)
