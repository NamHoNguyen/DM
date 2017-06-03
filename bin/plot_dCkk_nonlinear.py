import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3

fig = plt.figure(figsize=(10,8))
dash = [10,5]
suffix = ['_cut_ibarrier_iconc']
prefix = ['May21']
mass = [1.0]
label = ['Non-linear FDM']
for i in range(len(suffix)):
    ClsCDM = np.loadtxt('data/'+prefix[i]+'_matter2lens_WF_CDM'+suffix[i]+'_fCls.csv')
    Cls = np.loadtxt('data/'+prefix[i]+'_matter2lens_WF_FDM_'+str(mass[i])+suffix[i]+'_fCls.csv')
    plt.plot(Cls[:,0],(Cls[:,1]-ClsCDM[:,1])/ClsCDM[:,1],'navy',label=label[i])

##### LINEAR #####
ClFileCDM = 'data/Nov3_highAcc_CDM_unlensed_axion_fCls.csv'
ClFile1 = 'data/Nov3_highAcc_FDM_e-22_unlensed_axion_fCls.csv'
kkcol = 4

ClsCDM = np.loadtxt(ClFileCDM,delimiter=',')
Cls1 = np.loadtxt(ClFile1,delimiter=',')

lmin = 0
x = range(len(ClsCDM))
y1 = (Cls1[lmin:,kkcol]-ClsCDM[lmin:,kkcol])/ClsCDM[lmin:,kkcol]
l, = plt.plot(x[lmin:],y1,'navy',label='Linear FDM')
l.set_dashes(dash)

plt.xlabel('$L$',fontsize=24)
plt.ylabel('$(C^{\kappa\kappa}_L - C^{\kappa\kappa,\\rm CDM}_L)/ C^{\kappa\kappa,\\rm CDM}_L$',fontsize=24)
plt.xlim([0,50000])
plt.ylim([-1.,0.2])
plt.legend(loc='lower left')
plt.show()
#plt.savefig('output/June2_dCkk_nonlinear.pdf')    
