from driver import FisherForecast as clFish
import os
import numpy as np

iniFile = 'input/degeneracyTest.ini'
prefix = ''
dClsRoot ={}
nlkkLoc = {}
noiseTOverride={}
noise = [0.1]
beam = 0.3
#noise = [0.05,0.1,0.5]
#data = np.zeros([len(noise),2])
nlkkLoc['S4'] = 'data/Oct20_gradCut_2000_polComb_mv_beamY_2.0_noiseY_1.0_tellminY_30.0_tellmaxY_4000.0_kmax_5000.txt'
#nlkkLoc['S4'] = 'data/June13_gradCut_2000_polComb_mv_beamY_2.0_noiseY_0.1_tellminY_30.0_tellmaxY_4000.0_kmax_5000.txt'
#nlkkLoc['S4'] = 'data/June13_pmax10k_gradCut_2000_polComb_mv_beamY_2.0_noiseY_0.1_tellminY_30.0_tellmaxY_10000.0_kmax_10000.txt'
noiseTOverride['S4'] = 1.0
for i in range(len(noise)):
    source = 'CMB_noise'+str(noise[i])
    #nlkkLoc['CSST'] = 'data/June13_gradCut_2000_polComb_TT_beamY_'+str(beam)+'_noiseY_'+str(noise[i])+'_tellminY_100_tellmaxY_45000_kmax_60000.txt'
    #noiseTOverride['CSST'] = noise[i]
    nlkkLoc['DMS'] = 'data/June13_gradCut_2000_polComb_mv_beamY_0.167_noiseY_0.1_tellminY_30.0_tellmaxY_45000.0_kmax_60000.txt'
    F = clFish(iniFile,prefix=prefix,source=source,dClsRoot=dClsRoot,nlkkLocOverride=nlkkLoc,noiseTOverride=noiseTOverride)
    print "Calculating Fisher matrix..."
    FisherMat = F.calcFisher(verbose = True)
    #F.confEllipse('mnu','m_ax',confLevel=1,savefig=True,savedata=True,verbose=True)        
    #data[i,0] = noise[i]
    #data[i,1] = F.margSigma('m_ax')
    #data[i,1] = F.margSigma('mnu')
    #print data
#np.savetxt('sigma_m_ax_tau0.002_S4+CSST.csv',data)
#np.savetxt('sigma_mnu_tau0.002_S4+CSST+BAO',data)
                                                
