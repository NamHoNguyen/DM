import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
import itertools
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
dash1 = [10,5,2,5]
#dash = [[10,5],[10,5,2,5]]
color=itertools.cycle(['k','silver'])

fig = plt.figure(figsize=(9,7))

ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'
Cls1 = np.loadtxt(ClFile1)
l, = plt.plot(Cls1[:,0],Cls1[:,1],'m',label='$10^{-22}$eV FDM')

ClFile = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
Cls = np.loadtxt(ClFile)
l, = plt.plot(Cls[:,0],Cls[:,1],'c--',label='CDM')
l.set_dashes(dash1)

cdm = interp1d(Cls[:,0],Cls[:,1],bounds_error=True)
fdm = interp1d(Cls1[:,0],Cls1[:,1],bounds_error=True)

nbin = 20
Lrange = np.linspace(10000,60000,nbin)
#dataFile = ['June2_SNvsL_noise0.1.csv','June2_SNvsL_noise0.5.csv']
dataFile = ['June3_SNvsL_noise0.1.csv','June3_SNvsL_noise0.5.csv']
label = ['0.1 $\mu$K, 4000 deg$^2$','0.5 $\mu$K, 4000 deg$^2$']
ls = ['-','--']
shiftFactor = 0.
for index in range(len(dataFile)):
    data = np.loadtxt('data/'+dataFile[index])
    
    xs = []
    ys = []
    x = 0.
    y = 0.
    j = 1
    count = 0
    for i in range(len(data)):
        if Lrange[j] < data[i,0]:
            if count != 0:
                xs.append(x/count)
                ys.append(np.sqrt(y))
                count = 0
            j += 1
            x = 0.
            y = 0.
        count += 1
        x += data[i,0]
        y += data[i,1]**2
    if count != 0:
        xs.append(x/count)
        ys.append(np.sqrt(y))
        
    print 'Total SN: ',np.sqrt(np.sum(np.array(ys)**2))
    xs = np.array(xs)
    ys = np.array(ys)
    yerr = (cdm(xs)-fdm(xs))/ys
    (_,caps,eb)=plt.errorbar(xs+shiftFactor,fdm(xs+shiftFactor),yerr=yerr,ls='',ecolor=color.next(),capsize=5,label=label[index])
    for cap in caps:
        cap.set_markeredgewidth(2)
    eb[0].set_linestyle(ls[index])
    shiftFactor+=300

    
plt.xlabel('$L$',size=24)
plt.ylabel('$C_L^{\kappa\kappa}$',size=24)
plt.legend(loc='upper right',ncol=1)
plt.xlim([10000,45000])
plt.ylim([1.3e-11,2.2e-10])
plt.yscale('log')
#plt.show()
plt.savefig('output/June3_errorbar.pdf')
