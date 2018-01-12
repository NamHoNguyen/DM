from scipy.interpolate import interp1d
import numpy as np


for ClsFile in ['data/Aug6_highAcc_CDM_scalCls','data/Aug6_highAcc_CDM_lensedCls']:
    Cls = np.loadtxt(ClsFile+'.dat')
    ClsFile_new = ClsFile+'_new.dat'

    ell_new = np.arange(min(Cls[:,0]),max(Cls[:,0]),1)
    Cls_new = ell_new.reshape([len(ell_new),1])
    for i in range(1,len(Cls[0,:])):
        f = interp1d(Cls[:,0],Cls[:,i],kind='linear',bounds_error=False,fill_value=0.)
        Cls_new = np.hstack([Cls_new,f(ell_new).reshape([len(ell_new),1])])
    np.savetxt(ClsFile_new,Cls_new)
