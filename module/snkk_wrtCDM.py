import numpy as np
from sn_Mat import LensForecast

def snkk_wrtCDM(NlFile,ClsCDM,Cls1,Lmin,Lmax,dL,fsky):
    specType = 'kk'
    Nls = np.loadtxt(NlFile)
    # Change nan and inf to a number
    for Nli in range(len(Nls)):
        if (Nls[Nli,1] != Nls[Nli,1]) or (Nls[Nli,1] == float(str('inf'))):
            Nls[Nli,1] = 1.e300
            
    LFCDM = LensForecast()
    LF1 = LensForecast()
    
    LFCDM.loadKKDirect(ClsCDM[:,0],ClsCDM[:,1],Nls[:,0],Nls[:,1])
    LF1.loadKKDirect(Cls1[:,0],Cls1[:,1],Nls[:,0],Nls[:,1])
    #LF2.loadKKDirect(range(len(Cls2)),Cls2[:,kkcol],Nls[:,0],Nls[:,1])
    
    # L range
    ellBinEdges = np.arange(Lmin-dL/2.,Lmax+dL/2.,dL)
    ellMid  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
    #sigs = fsky*(2.*ellMid+1.)/2.*dL
    #print 'Predicted: < ',np.sqrt(sigs.sum())
    
    varCDM,sigsCDM,dummy = LFCDM.KnoxCov(specType,specType,ellBinEdges,fsky)
    var1,sigs1,dummy = LF1.KnoxCov(specType,specType,ellBinEdges,fsky)
    
    # Change denominator of S/N to CDM model
    
    ratio = var1/varCDM
    ratio[np.isnan(ratio)] = 1.
    sigs1 = sigs1*ratio
    
    # Get SN then take the rms
    SN1 = (np.sqrt(sigs1)-np.sqrt(sigsCDM))**2
    #SN1 = np.sqrt(SN1.sum())
    return ellMid,SN1
