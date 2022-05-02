# ----------------------------------------------------------------------
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
import h5py
#----------------------------------

# --- Plotting ---    

fig,ax = plt.subplots(5,1,sharex=True)    

ax[-1].set_xlabel('time',fontsize=18)

ax[0].set_ylabel(r'E(t)', fontsize=18)
ax[1].set_ylabel(r'$E$', fontsize=18)
ax[2].set_ylabel(r'$n_\alpha$', fontsize=18)
ax[3].set_ylabel(r'$J_{x,y}$', fontsize=18)
ax[4].set_ylabel(r'$D$', fontsize=18)

for i in range(0,len(ax)):
    ax[i].tick_params(labelsize=14,direction='in')
    ax[i].yaxis.set_ticks_position('both')
    ax[i].xaxis.set_ticks_position('both')

for i in range(1,len(sys.argv)):

    flin = sys.argv[i]
    f = h5py.File(flin, 'r')
 
    tp = np.array(f['time'])
    if f.attrs['applyfield'] == 1:
        EF = np.array(f['efield'])
        
    Ekin = np.array(f['ekin'])
    if 'etot' in f.keys():
        Etot = np.array(f['etot'])
    else:
        Etot = Ekin
    occ = f['occ']
    curr = f['current']
    if 'dipole' in f.keys():
        dip = np.array(f['dipole'])
        

    if f.attrs['applyfield'] == 1:
        ax[0].plot(tp,EF[:,0],c='purple',linewidth=2.0)
        ax[0].plot(tp,EF[:,1],c='green',linewidth=2.0)
  
    ax[1].plot(tp,Ekin,linewidth=2.0)
    ax[1].plot(tp[0:],Etot[0:],linewidth=2.0)

    for ibnd in range(occ.shape[1]):
        ax[2].plot(tp,occ[0:,ibnd],linewidth=2.0)

    ax[3].plot(tp,curr[0:],linewidth=2.0)

    if 'dipole' in f.keys():
        ax[4].plot(tp,dip[0:],linewidth=2.0)
    
    f.close()
    
plt.show()

exit()
