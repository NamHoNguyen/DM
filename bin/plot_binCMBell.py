import numpy as np
from snkk_wrtCDM import snkk_wrtCDM
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
dash=[10,5]

fig = plt.figure(figsize=(8,7))

ClFileCDM = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'

ClsCDM = np.loadtxt(ClFileCDM)
Cls1 = np.loadtxt(ClFile1)

dL = 100
beams = [0.3]
fskys = [0.1]
noises = [0.1,0.5]
Lmins = [100]
Lmaxs = [60000]
prefixs = ['Aug6_gradCut_2000_']
ests = ['TT']
kmaxs = [60000]
lRange = zip(np.arange(0,45000,500),np.arange(500,45500,500))
print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'

# Get Total S/N
lRange = [(100,45000)]
dataSum = np.zeros(len(noises))
i=0
for noise,(lmin,lmax),Lmax,beam,fsky,Lmin,prefix,est in itertools.product(noises,lRange,Lmaxs,beams,fskys,Lmins,prefixs,ests):
        NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
        ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
        SN1 = np.sqrt(SN1.sum())
        dataSum[i] = SN1
        print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
        i+=1
                                                    
##### Plot It #########
data1 = np.loadtxt('data/Aug6_binCMBell_beam0.3_noise0.1.csv')
data2 = np.loadtxt('data/Aug6_binCMBell_beam0.3_noise0.5.csv')
print 'RMS/Total=',[np.sqrt(sum(data1[:,1]**2)),np.sqrt(sum(data2[:,1]**2))]/dataSum
plt.plot(data1[:,0],data1[:,1]/dataSum[0]*100.,'k',label='0.1 $\mu$K-arcmin')
l, =plt.plot(data2[:,0],data2[:,1]/dataSum[1]*100.,'silver',ls='--',label='0.5 $\mu$K-arcmin')
l.set_dashes(dash)
plt.xlabel('$\ell$',fontsize=24)
plt.ylabel('% of Signal-to-Noise',fontsize=24)
plt.xlim([0,45000])
plt.legend(loc='upper right')
plt.locator_params(axis='x', nbins=5)
plt.savefig('output/Aug6_binCMBell_beam0.3.pdf')
#plt.show()
