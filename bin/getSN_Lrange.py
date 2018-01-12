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
prefixs = ['Aug6_gradCut_2000_']
ests = ['TT']
lmins = [100]
lmaxs= [45000]
kmaxs = [60000]
lRange = [[lmin,lmax] for lmin in lmins for lmax in lmaxs]
print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'
############## LMIN ####################
Lmins = np.linspace(1000,45000,30)
Lmaxs = [50000]
data = np.zeros([len(Lmins),2])
i = 0

for (lmin,lmax),Lmax,beam,fsky,noise,Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys,noises,Lmins,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'    
    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    SN1 = np.sqrt(SN1.sum())
    data[i,0] = Lmin
    data[i,1] = SN1
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1
np.savetxt('data/Aug6_SNvsLMin.csv',data)
############## LMAX ####################
Lmins = [100]
Lmaxs = np.linspace(1000,50000,30)
data = np.zeros([len(Lmaxs),2])
i = 0

for (lmin,lmax),Lmax,beam,fsky,noise,Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys,noises,Lmins,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    SN1 = np.sqrt(SN1.sum())
    data[i,0] =Lmax
    data[i,1] =SN1
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1
np.savetxt('data/Aug6_SNvsLMax.csv',data)
                        
