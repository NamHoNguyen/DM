import numpy as np
from scipy.interpolate import interp1d
from ConfigParser import SafeConfigParser
from orphics.tools.io import getFileNameString

def lensNoise(Config,expName,lensName,beamOverride=None,noiseTOverride=None,lkneeTOverride=None,lkneePOverride=None,alphaTOverride=None,alphaPOverride=None,tellminOverride=None,pellminOverride=None,tellmaxOverride=None,pellmaxOverride=None,kellmaxOverride=None):

    from orphics.tools.io import dictFromSection, listFromConfig

    beam = listFromConfig(Config,expName,'beams')
    noise = listFromConfig(Config,expName,'noises')
    freq = listFromConfig(Config,expName,'freqs')
    lkneeT,lkneeP = listFromConfig(Config,expName,'lknee')
    alphaT,alphaP = listFromConfig(Config,expName,'alpha')
    tellmin,tellmax = listFromConfig(Config,expName,'tellrange')
    if tellminOverride is not None: tellmin = tellminOverride
    if tellmaxOverride is not None: tellmax = tellmaxOverride
    pellmin,pellmax = listFromConfig(Config,expName,'pellrange')
    if pellminOverride is not None: pellmin = pellminOverride
    if pellmaxOverride is not None: pellmax = pellmaxOverride
    lmax = int(Config.getfloat(expName,'lmax'))

    pols = Config.get(lensName,'polList').split(',')
    freq_to_use = Config.getfloat(lensName,'freq')
    ind = np.where(np.isclose(freq,freq_to_use))
    beamFind = np.array(beam)[ind]
    noiseFind = np.array(noise)[ind]
    assert beamFind.size==1
    assert noiseFind.size==1
    if beamOverride is not None:
        beamX = beamY = beamOverride
    else:
        beamX = beamY = beamFind[0]
    if noiseTOverride is not None:
        noiseTX = noiseTY = noiseTOverride
    else:
        noiseTX = noiseTY = noiseFind[0]
    if lkneeTOverride is not None: lkneeT = lkneeTOverride
    if lkneePOverride is not None: lkneeP = lkneePOverride
    if alphaTOverride is not None: alphaT = alphaTOverride
    if alphaPOverride is not None: alphaP = alphaPOverride

    #import flipper.liteMap as lm
    import liteMap as lm
    from alhazen.quadraticEstimator import NlGenerator,getMax
    deg = 6. #for lmin=30
    px = 1.0
    dell = 10
    gradCut = 2000
    kellmin = 10
    lmap = lm.makeEmptyCEATemplate(raSizeDeg=deg, decSizeDeg=deg,pixScaleXarcmin=px,pixScaleYarcmin=px)
    kellmax = kellmaxOverride #max(tellmax,pellmax)
    from orphics.theory.cosmology import Cosmology
    cc = Cosmology(lmax=int(kellmax),pickling=True)
    theory = cc.theory
    bin_edges = np.arange(kellmin,kellmax,dell)
    myNls = NlGenerator(lmap,theory,bin_edges,gradCut=gradCut)
    myNls.updateNoise(beamX,noiseTX,np.sqrt(2.)*noiseTX,tellmin,tellmax,pellmin,pellmax,beamY=beamY,noiseTY=noiseTY,noisePY=np.sqrt(2.)*noiseTY,lkneesX=(lkneeT,lkneeP),lkneesY=(lkneeT,lkneeP),alphasX=(alphaT,alphaP),alphasY=(alphaT,alphaP))

    lsmv,Nlmv,ells,dclbb,efficiency = myNls.getNlIterative(pols,kellmin,kellmax,tellmax,pellmin,pellmax,dell=dell,halo=True)

    fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','tellminY','tellmaxY','kmax'],[gradCut,'mv',beamY,noiseTY,tellmin,tellmax,kellmax])+".txt"
    np.savetxt(fileName,np.vstack((lsmv,Nlmv)).transpose())
    return lsmv,Nlmv,ells,dclbb,efficiency

###########################################
iniFile = 'input/noise_S4.ini'
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)


saveRoot = "output/dump/May29" #"output/dump/Mar13"
#expName = 'DM-9.5arcsec' #'DM-18arcsec'
expName = 'S4'
lensName = 'lensing'
kmax = 6000

lmin = 100
lmax = 45000
#delensTolerance = 1.0

#tRange1 = [(x,y) for x in [100] for y in np.hstack([[1000],np.arange(5000,45000,1000)])]
#tRange2 = [(x,y) for x in np.hstack([[100,500],np.arange(1000,45000,1000)]) for y in [45000]]
#tRange = zip(np.arange(0,45000,500),np.arange(500,45500,500))
#tRange = np.vstack([tRange1,tRange2,tRange3])

#tRange = [(100,45000)]
lsmv,Nlmv,ells,dclbb,efficiency = lensNoise(Config,expName,lensName,kellmaxOverride=kmax)
#print lsmv,Nlmv,ells,dclbb,efficiency
