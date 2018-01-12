from driver import FisherForecast as clFish
import os
import numpy as np

iniFile = 'input/floorTest.ini'
prefix = ''
dClsRoot ={}
nlkkLoc = {}
noiseTOverride={}
noise = [0.05,0.1,0.5]
#noise = [2.5,2.0,1.0]
data = np.zeros([len(noise),2])
nlkkLoc['S4'] = 'data/Oct20_gradCut_2000_polComb_mv_beamY_2.0_noiseY_1.0_tellminY_30.0_tellmaxY_4000.0_kmax_5000.txt'
noiseTOverride['S4'] = 1.0
for i in range(len(noise)):
    nlkkLoc['CSST'] = 'data/Aug6_gradCut_2000_polComb_TT_beamY_0.3_noiseY_'+str(noise[i])+'_tellminY_100_tellmaxY_45000_kmax_60000.txt'
    noiseTOverride['CSST'] = noise[i]
    F = clFish(iniFile,prefix=prefix,dClsRoot=dClsRoot,nlkkLocOverride=nlkkLoc,noiseTOverride=noiseTOverride)
    print "Calculating Fisher matrix..."
    FisherMat = F.calcFisher(verbose = True)
    data[i,0] = noise[i]
    data[i,1] = F.margSigma('m_ax')
    #data[i,1] = F.margSigma('mnu')
    print data
np.savetxt('sigma_m_ax_tau0.002_S4+CSST.csv',data)
#np.savetxt('sigma_mnu_tau0.002_S4+CSST+BAO',data)
                                                
