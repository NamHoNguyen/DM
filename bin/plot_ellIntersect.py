import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'mathtext.fontset':'cm'})
#mpl.rcParams.update({'font.size': 18})

i = 0 # index for TT
Cls = np.loadtxt('../pyfisher/output/June7_newoptimal_vhhAcc_lensed_scalar_fCls.csv',delimiter=',')
ellCls = range(len(Cls[:,0]))
plt.loglog(ellCls,Cls[:,0])

for deproj in ['0','3']:
    NlTT = np.loadtxt('../pyfisher/tests/TT/SOV3_T_default1-4-2_noisecurves_deproj'+deproj+'_SENS1_mask_04000_ell_TT_yy.txt')
    plt.loglog(NlTT[:,0],NlTT[:,1],label='deproj'+deproj)

plt.axvline(x=2000,color='k')
plt.xlim([40,8000])
plt.ylim([1e-7,1e2])
plt.xlabel('$\ell$',fontsize=18)
plt.ylabel('$C_{\ell}^{TT}$',fontsize=18)
plt.legend()
plt.show()
