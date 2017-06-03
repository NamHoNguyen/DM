import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3

fig = plt.figure(figsize=(9,7))
size = 18
dash1 = [10,5,2,5]

delens_lmax = 2000
prefix = ['data/May27_gradCut_2000_polComb_']
suffix = ['_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_45000_kmax_60000','_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_39200_kmax_39200']
ests = ['TT','EE','TE','TB','EB'] #,'mv']

Nls = np.loadtxt(prefix[0]+'TT'+suffix[0]+'.txt')
plt.loglog(Nls[:,0],Nls[:,1],'k')
pos = 1400; a = 2.2
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'TT',size=size,color='k') #,weight='ultralight')

#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_EE_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000.txt')
Nls = np.loadtxt(prefix[0]+'EE'+suffix[0]+'.txt')
plt.loglog(Nls[:,0],Nls[:,1],color='g')
pos = 50; a = 1.5
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'EE',size=size,color='g') #,weight='ultralight')

#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TE_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000.txt')
Nls = np.loadtxt(prefix[0]+'ET'+suffix[0]+'.txt')
plt.loglog(Nls[:,0],Nls[:,1],'teal')
pos = 50; a = 1.5
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'ET',size=size,color='teal') #,weight='ultralight')

#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TB_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000.txt')
Nls = np.loadtxt(prefix[0]+'TB'+suffix[0]+'.txt')
plt.loglog(Nls[:,0],Nls[:,1],color='grey')
pos = 0; a = 1.1
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'TB',size=size,color='grey') #,weight='ultralight')
#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TB_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000_delens_1.0.txt')
Nls = np.loadtxt(prefix[0]+'TB'+suffix[1]+'_delens_1.0.txt')
plt.loglog(Nls[:delens_lmax,0],Nls[:delens_lmax,1],color='grey',ls='--')
pos = 1; a = 0.45
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'TB delensed',size=size,color='grey') #,weight='ultralight')

#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_EB_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000.txt')
Nls = np.loadtxt(prefix[0]+'EB'+suffix[0]+'.txt')
plt.loglog(Nls[:,0],Nls[:,1],'indianred')
pos = 3; a = 1.2
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'EB',size=size,color='indianred') #,weight='ultralight')
#Nls = np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_EB_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_30000_kmax_30000_delens_1.0.txt')
Nls = np.loadtxt(prefix[0]+'EB'+suffix[1]+'_delens_1.0.txt')
plt.loglog(Nls[:delens_lmax,0],Nls[:delens_lmax,1],'indianred',ls='--')
pos = 4; a = 1.3
plt.text(Nls[:,0][pos],a*Nls[:,1][pos],'EB delensed',size=size,color='indianred') #,weight='ultralight')

ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'
Cls1 = np.loadtxt(ClFile1)
#plt.loglog(range(len(Cls1)),Cls1[:,4],color='gray',ls='-.',label='FDM')
l, = plt.loglog(Cls1[:,0],Cls1[:,1],'m',label='$10^{-22}$eV FDM')

ClFile = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
Cls = np.loadtxt(ClFile)
l,= plt.loglog(Cls[:,0],Cls[:,1],'c--',label='CDM')
l.set_dashes(dash1)
pos = 30000; a = 0.15
plt.text(Cls[:,0][pos],a*Cls[:,1][pos],'$\kappa\kappa$',size=(size+8),color='m') #,weight='ultralight')

plt.xlim([30,1e5])
plt.ylim([5e-13,1e-5])
plt.xlabel('$L$',size=24)
plt.ylabel('$C_L^{\kappa\kappa}$ or $N_L^{\kappa\kappa}$',size=24)
plt.legend(loc='lower left',ncol=1)
plt.savefig('output/June2_noisecurves_beam0.3_noise0.1.pdf')
#plt.show()
