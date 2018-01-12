import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 10
dash = [[10,0.000001],[10,5],[10,5,2,5]]

fig = plt.figure(figsize=(8,7))

color = ['b','g','r']
#ls = ['-','--','-.']

exp = '_S4+CSST'
#exp = ''
data1 = np.flipud(np.loadtxt('data/sigma_m_ax_tau0.01'+exp+'.csv'))
data2 = np.flipud(np.loadtxt('data/sigma_m_ax_tau0.005'+exp+'.csv'))
data3 = np.flipud(np.loadtxt('data/sigma_m_ax_tau0.002'+exp+'.csv'))

noise = data1[:,0]
tau = [0.01,0.005,0.002]
for i in range(len(noise)):
    print 1.e-22/data1[i,1]
    sig_ma = np.array([data1[i,1],data2[i,1],data3[i,1]])*1.e22
    l, = plt.plot(tau,sig_ma,color[i]+'--o',label=str(noise[i])+' $\mu$K-arcmin')
    l.set_dashes(dash[i])    
plt.xlabel('$\sigma(\\tau)$',fontsize=19)
plt.ylabel('$\sigma(m_{\\rm FDM})$ [10$^{-22}$eV]',fontsize=19)
plt.legend(loc='center right')
plt.xlim([0.001,0.011])
plt.ylim([0.,0.15])
plt.savefig('output/Oct27_floorTest_mFDM'+exp+'.pdf')
#plt.show()
'''
exp = '_S4+CSST+BAO'
data1 = np.loadtxt('data/sigma_mnu_tau0.01'+exp)
data2 = np.loadtxt('data/sigma_mnu_tau0.005'+exp)
data3 = np.loadtxt('data/sigma_mnu_tau0.002'+exp)
data = [data1,data2,data3]
noise = data1[:,0]
tau = [0.01,0.005,0.002]
for i in range(len(data)):
    print 0.06/data[i][:,1]
    plt.plot(data[i][:,0],data[i][:,1]*1.e3,color[i]+ls[i]+'o',label='$\sigma(\\tau)='+str(tau[i])+'$')
plt.xlabel('noise')
plt.ylabel('$\sigma(m_{\\nu})$')
plt.legend(loc='upper right')
plt.xlim([0.04,0.51])
plt.ylim([0.,60.])
plt.savefig('output/Oct23_floorTest_mnu'+exp+'.png')
#plt.show()
'''
