import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
#.binned_statistic as binit
Cls = np.loadtxt('data/Aug6_highAcc_CDM_lensedCls.dat')[:,:2]
print Cls
i = 5000
Cls_raw = Cls[i:,:].copy()
plt.loglog(Cls_raw[:,0],Cls_raw[:,1])

dl = 50
nbin = int((max(Cls_raw[:,0])-min(Cls_raw[:,0]))/dl)
print nbin
#print Cls_raw[:,0]
val,ellBinEdges,dummy = stats.binned_statistic(Cls_raw[:,0],Cls_raw[:,1],statistic='mean',bins=nbin)
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
N = len(val)
Cls_raw = np.hstack([ellMids.reshape(N,1),val.reshape(N,1)])
plt.loglog(Cls_raw[:,0],Cls_raw[:,1])
plt.show()
