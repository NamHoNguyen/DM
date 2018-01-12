import numpy as np
from snkk_wrtCDM import snkk_wrtCDM
import matplotlib.pyplot as plt
import itertools
import matplotlib
matplotlib.rcParams.update({'font.family':'serif'})

ClFileCDM = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'

ClsCDM = np.loadtxt(ClFileCDM)
Cls1 = np.loadtxt(ClFile1)

#NlFile = 'data/BlakeNoisex2.csv'
#NlFile = 'data/Aug6_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_45000_kmax_60000.txt'
NlFile = 'dump/Sep22_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.5_tellminY_100_tellmaxY_45000_kmax_60000.txt'
dL = 1
Lmin = 100
Lmax = 60000
fsky = 0.1

ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
SN1 = np.sqrt(SN1.sum())
print SN1
#print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
