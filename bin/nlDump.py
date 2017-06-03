from orphics.tools.io import Plotter,dictFromSection,listFromConfig,getFileNameString
import liteMap as lm
#from szlib.szcounts import ClusterCosmology
#from alhazen.halos import NFWMatchedFilterSN
import numpy as np
from orphics.tools.cmb import loadTheorySpectraFromCAMB
from alhazen.quadraticEstimator import NlGenerator,getMax
import sys,os
import itertools

saveRoot = "output/dump/May30_nlDump" #"output/dump/Mar13"
mapRoot = os.environ['DERIVGEN_DIR']+"data/templateMap_27k.fits"
#mapRoot = os.environ['DERIVGEN_DIR']+"data/order5_lensedCMB_T_beam_cutout_3.fits"
# Make a CMB Noise Curve
#cambRoot = os.environ['CMBPROJ_DIR']+"data/ell28k_highacc"
cambRoot = os.environ['CMBPROJ_DIR']+"data/TheorySpectra/Nov10_highAcc_CDM"
#beamRange = [0.16] #TolTEC
beamRange = [2.0]
halo = True
#delensTolerance = None
delensTolerance = 1.0

noiseY = float(sys.argv[1])
#tellminY = int(sys.argv[2])
#pellminY = int(sys.argv[3])
tellminYRange = ([100])
#tellminYRange = ([500,1000,2000,3000,4000,5000])
#tellmaxYRange = ([1000,5000,10000,15000,20000,25000,30000])
#tellmaxYRange = ([30000])
tellmaxYRange = ([45000])

'''
tRange1 = [(x,y) for x in [100] for y in np.hstack([[1000],np.arange(5000,45000,1000)])]
tRange2 = [(x,y) for x in np.hstack([[100,500],np.arange(1000,45000,1000)]) for y in [45000]]
tRange3 = zip(np.arange(0,45000,500),np.arange(500,45500,500))
tRange = np.vstack([tRange1,tRange2,tRange3])
'''
#tRange = zip(np.arange(0,30000,2000),np.arange(2000,32000,2000))
#tRange = zip(np.arange(0,30000,15000),np.arange(15000,45000,15000))
#tRange = [(x,y) for x in np.hstack([[100,500,1000],np.arange(5000,16000,5000)]) for y in [30000]]
#tRange = [(x,y) for x in np.arange(16000,31000,1000) for y in [30000]]
#tRange = [(x,y) for x in [100] for y in np.arange(6000,30000,1000)]

#tRange = itertools.product(tellminYRange,tellmaxYRange)
tRange = [(30,4000)]
#beamRange = np.arange(0.5,5.0,0.5)

beamscale = lambda b: np.sqrt(8.*np.log(2.))*60.*180./np.pi/b

kmin = 40
deg = 5.
#px = 0.2
#px = 0.5
dell = 10
lpad = 60000
kmax = 6000
bin_edges = np.arange(kmin,kmax,dell)+dell

theory = loadTheorySpectraFromCAMB(cambRoot,unlensedEqualsLensed=True,useTotal=False,lpad=lpad)
#lmap = lm.makeEmptyCEATemplate(raSizeDeg=deg, decSizeDeg=deg,pixScaleXarcmin=px,pixScaleYarcmin=px)
#print lmap.data.shape
#lmap= lm.liteMapFromFits(mapRoot)
tmin,tmax=tRange[0]
px = np.round(180./max(kmax,tmin,tmax)*60.*0.95,2) #1.0
print 'pixel scale',px
lmap = lm.makeEmptyCEATemplate(raSizeDeg=deg, decSizeDeg=deg,pixScaleXarcmin=px,pixScaleYarcmin=px)   #print theory

i=0
for gradCut in [2000]:
#for gradCut in [2000]:
#for gradCut in [10000]:
    myNls = NlGenerator(lmap,theory,bin_edges,gradCut=gradCut)
    for polComb in ['EB','TB']:
    #for polComb in ['TE','EE','EB','ET','TB']:
    #for polComb in ['TT']:
        for beamY in beamRange:
            beamell = beamscale(beamY)
            #for tellminY,tellmaxY in [(100,30000)]:
            #for tellmaxY,pellmaxY in [(beamell,beamell)]:
            #for tellminY,tellmaxY in itertools.product(tellminYRange,tellmaxYRange):
            #for tellminY,tellmaxY in zip(np.arange(0,30000,500),np.arange(500,30500,500)):
            for tellminY,tellmaxY in tRange:
                pellminY = tellminY
                pellmaxY = tellmaxY
                #for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(noiseY,beamY,"sameGrad",100,2000,100,2000)]:
                #for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(noiseY,beamY,"sameGrad",100,51000,100,51000)]:
                #for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(noiseY,beamY,"sameGrad",100,30000,100,30000)]:
                for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(noiseY,beamY,"sameGrad",100,45000,100,45000)]:
                #for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(noiseY,beamY,"sameGrad",tellminY,tellmaxY,pellminY,pellmaxY)]:
                #for noiseX,beamX,lab,tellminX,tellmaxX,pellminX,pellmaxX in [(0.5,0.3,"sameGrad",500,5000,500,5000)]:

                    noiseTX = noiseX
                    noisePX = np.sqrt(2.)*noiseTX
                    noiseTY = noiseY
                    noisePY = np.sqrt(2.)*noiseTY

                    #kmax = getMax(polComb,tellmaxY,pellmaxY)
                    #kmax = 51000
                    print kmax,polComb
                    i+=1
                    print i,tellminY,tellmaxY,kmax,"delens:",delensTolerance

                    myNls.updateBins(bin_edges)
                    myNls.updateNoise(beamX,noiseTX,noisePX,tellminX,tellmaxX, \
                                      pellminX,pellmaxX,beamY=beamY,noiseTY=noiseTY, \
                                      noisePY=noisePY,tellminY=tellminY,tellmaxY=tellmaxY, \
                                      pellminY=pellminY,pellmaxY=pellmaxY)



                    
                    ls,Nls = myNls.getNl(polComb=polComb,halo=halo)

                    fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','tellminY','tellmaxY','kmax'],[gradCut,polComb,beamY,noiseY,tellminY,tellmaxY,kmax])+".txt"
                    np.savetxt(fileName,np.vstack((ls,Nls)).transpose())

                    if (polComb=='EB' or polComb=='TB') and (delensTolerance is not None):
                        ls, Nls, efficiency = myNls.iterativeDelens(polComb,delensTolerance,halo)
                        print 'efficiency ', efficiency
                        fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','tellminY','tellmaxY','kmax','delens'],[gradCut,polComb,beamY,noiseY,tellminY,tellmaxY,kmax,delensTolerance])+".txt"
                        
                        #fileName = saveRoot + getFileNameString(['gradCut','polComb','beamY','noiseY','grad','tellminY','pellminY','tellmaxY','pellmaxY','kmin','deg','px','delens'],[gradCut,polComb,beamY,noiseY,lab,tellminY,pellminY,tellmaxY,pellmaxY,kmin,deg,px,delensTolerance])+".txt"
                        np.savetxt(fileName,np.vstack((ls,Nls)).transpose())
