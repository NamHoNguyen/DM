import numpy as np
from sn_forecast import LensingForecast
import itertools

ClsCDM = np.loadtxt('data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv')
ClsFDM = np.loadtxt('data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv')
LF = LensingForecast(ClsCDM,ClsFDM)

# L range
Lmin,Lmax,dL = 100,60000,100
ellBinEdges = np.arange(Lmin-dL/2.,Lmax+dL/2.,dL)
    
fskys = [0.1,0.025]
noises = [0.1,0.5]
beams = [0.3,0.16]
prefixs = ['Aug6_gradCut_2000_']
ests = ['TT']
lmins = [100]
lmaxs = [45000]
lRange = [[lmin,lmax] for lmin in lmins for lmax in lmaxs]

print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'
for (lmin,lmax),beam,fsky,noise,prefix,est in itertools.product(lRange,beams,fskys,noises,prefixs,ests):
    NlFile = 'data/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
    Nls = np.loadtxt(NlFile)
    sn2 = LF.SNR2_fromNl(Nls,ellBinEdges,fsky)
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,np.sqrt(np.sum(sn2))
