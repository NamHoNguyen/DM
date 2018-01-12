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

dL = 100
#specType = 'kk'
fskys = [0.1,0.05,0.025,0.0125,0.00625,0.004]
noises = [0.5,0.354,0.25,0.177,0.125,0.1]
fskys_noises = zip(fskys,noises)
beams = [0.3]
Lmins = [100]
Lmaxs = [60000]
#prefixs = ['May26_gradCut_2000_'] #,'May26_new_gradCut_10000_','May26_new_gradCut_45000_'] #,'Feb9_gradCut_10000_']
prefixs = ['Sep22_gradCut_2000_'] #,'May26_new_gradCut_10000_','May26_new_gradCut_45000_'] #,'Feb9_gradCut_10000_']
ests = ['TT']
lmins = [100]
lmaxs = [45000]
lRange = [[lmin,lmax] for lmin in lmins for lmax in lmaxs]

print 'lmin,lmax,Lmin,Lmax,beam,fksy,noise,prefix,est,S/N'
i = 0
for (lmin,lmax),Lmax,beam,(fsky,noise),Lmin,prefix,est in itertools.product(lRange,Lmaxs,beams,fskys_noises,Lmins,prefixs,ests):
    NlFile = 'dump/'+prefix+'polComb_'+est+'_beamY_'+str(beam)+'_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
    #NlFile = 'output/dump/'+prefix+'polComb_'+est+'_beamY_0.16_noiseY_'+str(noise)+'_tellminY_'+str(lmin)+'_tellmaxY_'+str(lmax)+'_kmax_60000.txt'
    #NlFile = os.environ['CMBPROJ_DIR']+'data/NoiseCurvesKK/S4_Lens_temp_mv.csv'

    ellMid,SN1 = snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky)
    SN1 = np.sqrt(SN1.sum())
    #print 'Factor of : ',0.1/fsky,(0.5/noise)**2
    print lmin,lmax,Lmin,Lmax,beam,fsky,noise,prefix,est,SN1 #,SN2,SN3
    i+=1
