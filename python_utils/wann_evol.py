import sys
import os
import numpy as np
import f90nml
#--------------------------------------------------------------------------------------
class WannierEvolution():
    #========================================
    def __init__(self,PathExe,PathInp='./inp/',PathOut='./out/',PathLog='./log/',mpicmd=""):
        if len(mpicmd) > 0:
            self.exe = PathExe + "/exe/wann_evol_mpi.x"
        else:    
            self.exe = PathExe + "/exe/wann_evol.x"
        self.mpi = True
        self.mpicmd = mpicmd
        self.PathInp = PathInp
        self.PathOut = PathOut
        self.PathLog = PathLog
    #========================================
    def SetSysParams(self,MuChem,Beta,file_ham,gauge,FixMuChem=True,Filling=1.0,Output_Occ_KPTS=False):
        self.sysparams = {
            'MuChem': MuChem,
            'Beta': Beta,
            'file_ham': file_ham,
            'gauge': gauge,
            'FixMuChem': FixMuChem,
            'Filling': Filling,
            'Output_Occ_KPTS': Output_Occ_KPTS
        }
    #========================================
    def SetTimeParams(self,Nt,Tmax,file_field="",output_step=1,tstart=0.0,\
        relaxation_dynamics=None,T1_relax=None,T2_relax=None):
        self.timeparams = {
            'Nt': Nt,
            'Tmax': Tmax,
            'file_field': file_field,
            'output_step': output_step,
            'tstart': tstart
        }

        if relaxation_dynamics != None:
            self.timeparams['relaxation_dynamics'] = relaxation_dynamics

        if T1_relax != None:
            self.timeparams['T1_relax'] = T1_relax
        if T2_relax != None:
            self.timeparams['T2_relax'] = T2_relax
    #========================================
    def SetKPTS(self,kpoints_type,file_kpts="",nk1=1,nk2=1,nk3=1):
        self.kpoints = {
            'kpoints_type': kpoints_type,
            'file_kpts': file_kpts,
            'nk1': nk1,
            'nk2': nk2,
            'nk3': nk3
        }
    #========================================
    def WriteInput(self,file_inp):
        inp = {
            'SYSPARAMS': self.sysparams,
            'TIMEPARAMS': self.timeparams,
            'KPOINTS': self.kpoints
            }
        with open(file_inp, 'w') as nml_file:
            f90nml.write(inp, nml_file)
    #========================================
    def Run(self,prefix,debug_mode=False):
        nt = self.timeparams['Nt']
        tmax = self.timeparams['Tmax']
        file_pref = prefix + "_nt{}_tmax{}".format(nt,tmax)
  
        file_inp = self.PathInp + file_pref + '.inp'
        file_out = self.PathOut + file_pref
        file_log = self.PathLog + file_pref + '.log'

        self.WriteInput(file_inp)

        if(self.mpi):
            if debug_mode:
                os.system(self.mpicmd + ' ' + self.exe + ' ' + file_inp + ' ' + file_out)
            else:
                os.system(self.mpicmd + ' ' + self.exe + ' ' + file_inp + ' ' + file_out + ' > ' + file_log)
        else:
            if debug_mode:
                os.system(self.exe + ' ' + file_inp + ' ' + file_out)
            else:
                os.system(self.exe + ' ' + file_inp + ' ' + file_out + ' > ' + file_log)
#--------------------------------------------------------------------------------------

