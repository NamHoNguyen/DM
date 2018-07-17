import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3

dash2 = [10,5]
dash = [[1,0.001],[2,5],[10,5,2,5]]

fig = plt.figure(figsize=(10,8))
c = ['navy','m','lime']
suffix = ['_cut_ibarrier_iconc','_cut_ibarrier_iconc','_cut_ibarrier_iconc']
#prefix = ['May21','May21','May21']
prefix = ['Apr3','Apr3','Apr3']
mass = ['_FDM_1.0','_FDM_0.01','_WDM_1.0']
label = ['$10^{-22}$eV FDM','$10^{-24}$eV FDM','1keV WDM']
for i in range(len(suffix)):
    print 'Working on '+label[i]
    ClsCDM = np.loadtxt('data/'+prefix[i]+'_matter2lens_WF_CDM'+suffix[i]+'_fCls.csv')
    Cls = np.loadtxt('data/'+prefix[i]+'_matter2lens_WF'+mass[i]+suffix[i]+'_fCls.csv')
    y = (Cls[:,1]-ClsCDM[:,1])/ClsCDM[:,1]
    l, = plt.plot(Cls[:,0],y,c[i],ls='--',label=label[i])
    l.set_dashes(dash[i])
##### BARYON #####
ClsDMONLY = np.loadtxt('data/May22_matter2lens_DMONLY_WMAP7_fCls.csv')
ClsAGN = np.loadtxt('data/May22_matter2lens_AGN_WMAP7_fCls.csv')
y = (ClsAGN[:,1]-ClsDMONLY[:,1])/ClsDMONLY[:,1]
l, = plt.plot(ClsAGN[:,0],y,'maroon',ls='--',label='CDM+Baryons')
l.set_dashes(dash2)

plt.xscale('log')
plt.legend(loc = 'lower left')
plt.xlabel('$L$',fontsize=24)
plt.ylabel('$(C^{\kappa\kappa}_L - C^{\kappa\kappa,\\rm CDM}_L)/ C^{\kappa\kappa,\\rm CDM}_L$',fontsize=24)
plt.xlim([10,1e8])
plt.show()
#plt.savefig('output/June2_dCkk_WF.pdf')
