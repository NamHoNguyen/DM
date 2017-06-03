import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
#import matplotlib.pyplot.set_dashes
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['lines.linewidth'] = 3
size = 18
dash1 = [10,5,2,5]

#ClFileCDM = os.environ['FISHER_DIR']+'output/Nov3_highAcc_CDM_unlensed_axion_fCls.csv'
#ClFile1 = os.environ['FISHER_DIR']+'output/Nov3_highAcc_WDM_e-22_unlensed_axion_fCls.csv'
#kkcol = 4

#ClsCDM = np.loadtxt(ClFileCDM,delimiter=',')
#Cls1 = np.loadtxt(ClFile1,delimiter=',')

lmin = 0
#print y2

Nls=np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.5_tellminY_100_tellmaxY_30000_kmax_30000.txt')
plt.loglog(Nls[:,0],Nls[:,1],'k-',label='0.5 $\mu$K-arcmin')
pos = 50; a = 1.2
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'TT',size=size,color='k') #,weight='ultralight') 
a=np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000.txt')
l, = plt.loglog(a[:,0],a[:,1],'k--',label='0.1 $\mu$K-arcmin')
#l.set_dashes(dash1)
#x = range(len(ClsCDM))
'''
y1 = (ClsCDM[lmin:,kkcol]-Cls1[lmin:,kkcol])
#y1 = 1.e-7*(1+(Cls1[lmin:,kkcol]-ClsCDM[lmin:,kkcol])/ClsCDM[lmin:,kkcol])
#y2 = 1.e-7*(1.+0./ClsCDM[lmin:,kkcol])
#plt.loglog(x,y1,'k',label='$\Delta C^{\kappa\kappa}_L/ C^{\kappa\kappa}_L$')
plt.loglog(x,y1,'k',label='$|\Delta C^{\kappa\kappa}_L|$')
x = np.linspace(1.e2,1.e5,100)
y2 = 0.*x + (ClsCDM[lmin:,kkcol]-Cls1[lmin:,kkcol])[len(ClsCDM[lmin:,kkcol])-1]
plt.loglog(x,y2,'k--')
plt.text(5000,max(y2)*0.4,'min($|\Delta C^{\kappa\kappa}_L|$)')
'''
ClFile = os.environ['FISHER_DIR']+'output/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
Cls = np.loadtxt(ClFile)
plt.loglog(Cls[:,0],Cls[:,1],'c',label='CDM')
pos = 35; a = 1.2
plt.text(Cls[:,0][pos],a*Cls[:,1][pos],'$\kappa\kappa$',size=(size+8),color='m') #,weight='ultralight')


#color,ls = cs.next()
#ClFile1 = os.environ['FISHER_DIR']+'output/Nov3_highAcc_WDM_e-22_unlensed_axion_fCls.csv'
ClFile1 = os.environ['FISHER_DIR']+'output/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'
Cls1 = np.loadtxt(ClFile1)
#plt.loglog(range(len(Cls1)),Cls1[:,4],color='gray',ls='-.',label='FDM')
l, = plt.loglog(Cls1[:,0],Cls1[:,1],color='m',ls='--',label='FDM, $m=10^{-22}$eV')
l.set_dashes(dash1)

plt.xlabel('$L$',size=20)
plt.ylabel('$N_L^{\kappa\kappa}$',size=20)
plt.xlim([30,1.e5])
plt.ylim([5e-13,1e-5])
#plt.xscale('log')
plt.legend(loc='lower left')
plt.savefig('output/May22_noiseCompare.pdf')
#plt.show()
