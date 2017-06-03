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
noises = [0.1,0.5]
Lmins = [100]
Lmaxs = [60000]
prefixs = ['May26_gradCut_2000_']
ests = ['TT']
lmins = [100]
lmaxs= [45000]
kmaxs = [60000]
lRange = [[lmin,lmax] for lmin in lmins for lmax in lmaxs]
print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'

i = 0

for (lmin,lmax),Lmax,beam,fsky,noise,Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys,noises,Lmins,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'    

    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    data = np.transpose(np.vstack([ellMid,np.sqrt(SN1)]))
    np.savetxt('data/June3_SNvsL_noise'+str(noise)+'.csv',data)
    SN1 = np.sqrt(SN1.sum())
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1


