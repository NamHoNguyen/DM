import time
import numpy as np
from ConfigParser import SafeConfigParser
from scipy.interpolate import interp1d
from orphics.tools.io import getFileNameString
from orphics.tools.cmb import loadTheorySpectraFromCAMB
import os

def lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=None,pellminOverride=None,tellmaxOverride=None,pellmaxOverride=None,kellmaxOverride=None,tRange=None,delensTolerance = None):

    from orphics.tools.io import dictFromSection, listFromConfig

    beam = listFromConfig(Config,expName,'beams')
    noise = listFromConfig(Config,expName,'noises')
    freq = listFromConfig(Config,expName,'freqs')
    lkneeT,lkneeP = listFromConfig(Config,expName,'lknee')
    alphaT,alphaP = listFromConfig(Config,expName,'alpha')
    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    if tellminOverride is not None: tellminY = tellminOverride
    else: tellminY = tellmin
    if tellmaxOverride is not None: tellmaxY = tellmaxOverride
    else: tellmaxY = tellmax
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')
    if pellminOverride is not None: pellminY = pellminOverride
    else: pellminY = pellmin
    if pellmaxOverride is not None: pellmaxY = pellmaxOverride
    else: pellmaxY = pellmax
    
    if tRange is None:
        tRange = [(tellminY,tellmaxY)]
    
    pols = Config.get(lensName,'polList').split(',')
    '''
    freq_to_use = Config.getfloat(lensName,'freq')
    ind = np.where(np.isclose(freq,freq_to_use))
    beamFind = np.array(beam)[ind]
    noiseFind = np.array(noise)[ind]
    assert beamFind.size==1
    assert noiseFind.size==1
    '''
    if beamOverride is not None:
        beamX = beamY = beamOverride
    else:
        beamRange = np.array(beam)
        #beamX = -1
        #beamX = beamY = beamFind[0]
    if noiseTOverride is not None:
        noiseTX = noiseTY = noiseTOverride
    else:
        noiseRange = np.array(noise)
        #noiseTX = -1
        #noiseTX = noiseTY = noiseFind[0]
    if lkneeTOverride is not None: lkneeT = lkneeTOverride
    if lkneePOverride is not None: lkneeP = lkneePOverride
    if alphaTOverride is not None: alphaT = alphaTOverride
    if alphaPOverride is not None: alphaP = alphaPOverride

    import liteMap as lm
    from alhazen.quadraticEstimator import NlGenerator,getMax
    deg = 5.
    dell = 10
    gradCuts = [2000] #10000
    kmin = 40
    #kmax = max(tellmax,pellmax)

    if kellmaxOverride is not None: kmax = kellmaxOverride
    
    px = np.round(180./max(kmax,tellmax,tellmaxY)*60.*0.95,2) #1.0
    print 'pixel scale',px
    lmap = lm.makeEmptyCEATemplate(raSizeDeg=deg, decSizeDeg=deg,pixScaleXarcmin=px,pixScaleYarcmin=px)
    #lmap = lm.liteMapFromFits(os.environ['DERIVGEN_DIR']+"data/templateMap_27k.fits")
    lpad = 60000
    
    #cambRoot = os.environ['CMBPROJ_DIR']+"data/TheorySpectra/Nov10_highAcc_CDM"
    cambRoot = "data/Aug6_highAcc_CDM"
    theory = loadTheorySpectraFromCAMB(cambRoot,unlensedEqualsLensed=True,useTotal=False,lpad=lpad)
    '''
    from orphics.theory.cosmology import Cosmology
    cc = Cosmology(lmax=int(kmax),pickling=True)
    theory = cc.theory
    '''          
    i = 0
    bin_edges = np.arange(kmin,kmax,dell)
    print gradCuts,pols,beamRange,noiseRange,tRange
    for gradCut in gradCuts:
        myNls = NlGenerator(lmap,theory,bin_edges,gradCut=gradCut)
        # remove hardcoded 9000 limit
        myNls.N.bigell = lpad
        myNls.N.lmax_T = lpad
        myNls.N.lmax_P = lpad
        
        for polComb in pols:
            for beamX in beamRange:
                beamY = beamX
                for noiseTX in noiseRange:
                    noiseTY = noiseTX
                    for tellminY,tellmaxY in tRange:
                        #pellmin = pellminY = tellminY
                        #pellmax = pellmaxY = tellmaxY
                
                        i+=1
                        print i,tellminY,tellmaxY,kmax,polComb,beamX,noiseTX,"delens:",delensTolerance
                        #myNls.updateBins(bin_edges)
                        myNls.updateNoise(beamX,noiseTX,np.sqrt(2.)*noiseTX,tellmin,tellmax,pellmin,pellmax,beamY,noiseTY,np.sqrt(2.)*noiseTY,tellminY,tellmaxY,pellminY,pellmaxY,lkneesX=(lkneeT,lkneeP),lkneesY=(lkneeT,lkneeP),alphasX=(alphaT,alphaP),alphasY=(alphaT,alphaP))
                        
                        #lsmv,Nlmv,ells,dclbb,efficiency = myNls.getNlIterative(pols,kellmin,kellmax,tellmax,pellmin,pellmax,dell=dell,halo=True)
                        start = time.time()
                        ls,Nls = myNls.getNl(polComb=polComb,halo=True)
                        end = time.time()
                        print 'getNl:',end-start,'sec'
                        
                        noiseY = noiseTY
                        fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','tellminY','tellmaxY','kmax'],[gradCut,polComb,beamY,noiseY,tellminY,tellmaxY,kmax])+".txt"
                        np.savetxt(fileName,np.vstack((ls,Nls)).transpose())
                        
                        if (polComb=='EB' or polComb=='TB') and (delensTolerance is not None):
                            ls, Nls, efficiency = myNls.iterativeDelens(polComb,delensTolerance,True)
                            print 'efficiency ', efficiency
                            fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','tellminY','tellmaxY','kmax','delens'],[gradCut,polComb,beamY,noiseY,tellminY,tellmaxY,kmax,delensTolerance])+".txt"
                            np.savetxt(fileName,np.vstack((ls,Nls)).transpose())
                    
    return 0
    #return lsmv,Nlmv,ells,dclbb,efficiency

###########################################
iniFile = 'input/new_nlDump.ini'
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)


saveRoot = "dump/Oct20" #"output/dump/Mar13"
#expName = 'DM-9.5arcsec' #'DM-18arcsec'
#expName = 'DM-18arcsec'
expName = 'DM-S4'
lensName = 'lensing'
kmax = 5000 
#kmax = 60000 
#kmax = 39200

#lmin = 100
#lmax = 45000
delensTolerance = None
#delensTolerance = 1.0

#tRange1 = [(x,y) for x in [100] for y in np.hstack([[1000],np.arange(5000,45000,1000)])]
#tRange2 = [(x,y) for x in np.hstack([[100,500],np.arange(1000,45000,1000)]) for y in [45000]]
#tRange3 = zip(np.arange(0,45000,500),np.arange(500,45500,500))
#tRange = np.vstack([tRange1,tRange2,tRange3])

#tRange = [(100,39200)]
#tRange = [(100,45000)]
#lensNoise(Config,expName,lensName,kellmaxOverride=kmax,tRange=tRange,delensTolerance = delensTolerance)
lensNoise(Config,expName,lensName,kellmaxOverride=kmax,delensTolerance = delensTolerance)
#lensNoise(Config,expName,lensName,tellminOverride=lmin,pellminOverride=lmin,tellmaxOverride=lmax,pellmaxOverride=lmax,kellmaxOverride=kmax,tRange=tRange,delensTolerance = delensTolerance)
#print 'Results: ',lsmv,Nlmv,ells,dclbb,efficiency
