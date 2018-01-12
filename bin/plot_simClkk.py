import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10

fig = plt.figure(figsize=(10,8))
dash = [10,5]

name = ['clkk_theory','inputXrecon','reconXrecon_raw','totPower','varPower']
label = ['','Input x Reconstruction','Total Autospectrum','Theory Total Power','Variance of Power']
ls = ['-','o','o','--','^']
c = ['k','#3498DB','#F0B27A','r','#28B463']
#c = ['k','b','#orange','r','g']
for i in range(len(name)):
    data = np.loadtxt('data/_plot_'+name[i]+'.txt')
    plt.semilogy(data[:,0],data[:,1],ls[i],color=c[i],label=label[i],markeredgewidth=0.)
                

#l, = plt.plot(x[lmin:],y1,'navy',label='Linear FDM')
#l.set_dashes(dash)

plt.xlabel('$L$',fontsize=24)
plt.ylabel('$C^{\kappa\kappa}_L$',fontsize=24)
plt.xlim([1000,39000])
plt.ylim([1.e-12,1.e-7])
plt.legend(loc='lower left')
#plt.show()
plt.savefig('output/Oct3_simClkk.png')    
