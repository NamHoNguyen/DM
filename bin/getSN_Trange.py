import numpy as np
from snkk_wrtCDM import snkk_wrtCDM
import itertools
import matplotlib as mpl

ClFileCDM = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'

ClsCDM = np.loadtxt(ClFileCDM)
Cls1 = np.loadtxt(ClFile1)

dL = 100
beams = [0.3]
fskys = [0.1]
noises = [0.1]
Lmins = [100]
Lmaxs = [60000]
prefixs = ['Aug6_gradCut_2000_']
ests = ['TT']
print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'
############## TMIN ####################
lRange = [(x,y) for x in np.hstack([[100,500],np.arange(1000,45000,1000)]) for y in [45000]]
data = np.zeros([len(lRange),2])
i = 0

for (lmin,lmax),Lmax,beam,fsky,noise,Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys,noises,Lmins,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'    
    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    SN1 = np.sqrt(SN1.sum())
    data[i,0] = lmin
    data[i,1] = SN1
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1
np.savetxt('data/Aug6_SNvsTMin.csv',data)
############## LMAX ####################
lRange = [(x,y) for x in [100] for y in np.hstack([[1000],np.arange(5000,46000,1000)])]
data = np.zeros([len(lRange),2])
i = 0

for (lmin,lmax),Lmax,beam,fsky,noise,Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys,noises,Lmins,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    SN1 = np.sqrt(SN1.sum())
    data[i,0] =lmax
    data[i,1] =SN1
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1
np.savetxt('data/Aug6_SNvsTMax.csv',data)
                        
