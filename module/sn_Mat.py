import matplotlib
matplotlib.use('Agg')
import numpy as np
import csv
from scipy.interpolate import interp1d
from math import pi
import sys

from cmbUtils import smartCls
from mmUtils import Plotter
from lensingVar import TheorySpectra
from ConfigParser import ConfigParser

        
class LensForecast:

    def __init__(self):
        '''
        Make S/N projections for CMB and OWL auto and cross-correlations.

        K refers to the CMB (source) kappa
        S refers to the shear/kappa of an optical background galaxy sample
        G refers to the number density of an optical foreground galaxy sample

        '''
        self._haveKK = False
        self._haveKG = False
        self._haveGG = False
        
        self._haveSS = False
        self._haveSG = False
        
        self._haveKS = False

        self.theory = TheorySpectra()
        self.Nls = {}
        
    def loadKK(self,filename,colnum,Nfilename,Ncolnum,norm="lsq",transpower=[4.,0.25],nnorm="none",ntranspower=[0.,1.0]):#,nnorm="lsq",ntranspower=[2.0,0.25]):
        '''
        Loads CMB kappa Cls and noise curve. Note that this function uses cmbutils.smartCls which
        determines if colnum = 0 contains ells and uses it, else assumes a range of 2 to ellmax.
        
        
        '''

        if colnum==None:
            colnum = 11
            print "WARNING: Using default column no. 11 for Clkk"
        if Ncolnum==None:
            Ncolnum = 1
            print "WARNING: Using default column no. 1 for Nlkk"

        
        
        
        readCls = smartCls(filename)
        Cls = np.array(readCls.getCol(colnum=colnum,norm=norm,transpower=transpower))
        readNls = smartCls(Nfilename)
        Nls = np.array(readNls.getCol(colnum=Ncolnum,norm=nnorm,transpower=ntranspower))

        
        self.Nls['kk'] = interp1d(readNls.ells,Nls,bounds_error=False,fill_value=np.inf)
        self.theory.loadGenericCls(readCls.ells,Cls,'kk')
    
        self._haveKK = True

    def loadKKDirect(self,ellsCls,Cls,ellsNls,Nls):
        self.Nls['kk'] = interp1d(ellsNls,Nls,bounds_error=False,fill_value=np.inf)
        self.theory.loadGenericCls(ellsCls,Cls,'kk')
    
        self._haveKK = True
        
    def loadGG(self,filename,colnum,ngal,norm="lsq"):

    
        if colnum==None:
            colnum = 16
            print "WARNING: Using default column no. 16 for Clgg"

        self.ngalForeground = ngal
        
        readCls = smartCls(filename)
        Cls = np.array(readCls.getCol(colnum=colnum,norm=norm))

        self.theory.loadGenericCls(readCls.ells,Cls,'gg')

        self.Nls['gg'] = lambda x: 1./(self.ngalForeground*1.18e7)

        self._haveGG = True

    def loadGGDirect(self,ellsCls,Cls,ngal):
        self.ngalForeground = ngal
        self.Nls['gg'] = lambda x: 1./(self.ngalForeground*1.18e7)
        self.theory.loadGenericCls(ells,ellsCls,'gg')
    
        self._haveGG = True
        
        
    def loadSS(self,filename,colnum,ngal,shapeNoise=0.3,norm="lsq"):

        # fix
        if colnum==None:
            colnum = 16
            print "WARNING: Using default column no. 16 for Clgg"



        if shapeNoise==None or shapeNoise<1.e-9:
            print "No/negligible shape noise given. Using default = 0.3."
            self.shapeNoise=0.3

        else:             
            self.shapeNoise = shapeNoise
        self.ngalBackground = ngal
        self.Nls['ss'] = lambda x: self.shapeNoise*self.shapeNoise/(2.*self.ngalBackground*1.18e7)

        readCls = smartCls(filename)
        Cls = np.array(readCls.getCol(colnum=colnum,norm=norm))

        self.theory.loadGenericCls(readCls.ells,Cls,'ss')
        
        self._haveSS = True

    def loadSSDirect(self,ellsCls,Cls,ngal,shapeNoise=0.3):
        if shapeNoise==None or shapeNoise<1.e-9:
            print "No/negligible shape noise given. Using default = 0.3."
            self.shapeNoise=0.3

        else:             
            self.shapeNoise = shapeNoise
        self.ngalBackground = ngal
        self.Nls['ss'] = lambda x: self.shapeNoise*self.shapeNoise/(2.*self.ngalBackground*1.18e7)


        self.theory.loadGenericCls(ellsCls,Cls,'ss')
        
        self._haveSS = True

    def loadSG(self,filename,colnums,norm="lsq"):

    
        if colnums==None:
            colnums = [21,25]
            print "WARNING: Using default column nos. 21 and 25 for Clsg"



        
        readCls = smartCls(filename)
        Cls1 = np.array(readCls.getCol(colnum=colnums[0],norm=norm))
        Cls2 = np.array(readCls.getCol(colnum=colnums[1],norm=norm))

        try:
            assert np.allclose(Cls1,Cls2)
        except AssertionError:
            print "Cls1s2 not the same as Cls2s1!"
            print Cls1
            print Cls2
            sys.exit(1)

        self.theory.loadGenericCls(readCls.ells,Cls1,'sg')
        
        self._haveSG = True
        

    def loadSGDirect(self,ellsCls,Cls):
        self.theory.loadGenericCls(ellsCls,Cls,'sg')
        
        self._haveSG = True

    def loadKG(self,filename,colnum,norm="lsq",transpower=[2.,0.5]):

        if colnum==None:
            colnum = 12
            print "WARNING: Using default column no. 12 for Clkg"

        
        
        
        readCls = smartCls(filename)
        Cls = np.array(readCls.getCol(colnum=colnum,norm=norm,transpower=transpower))
        self.theory.loadGenericCls(readCls.ells,Cls,'kg')
        self._haveKG = True

    def loadKGDirect(self,ellsCls,Cls):
        self.theory.loadGenericCls(ellsCls,Cls,'kg')
        self._haveKG = True
                
    def loadKS(self,filename,colnum,norm="lsq",transpower=[2.,0.5]):

        if colnum==None:
            colnum = 12
            print "WARNING: Using default column no. 12 for Clks"

        
        
        
        readCls = smartCls(filename)
        Cls = np.array(readCls.getCol(colnum=colnum,norm=norm,transpower=transpower))
        self.theory.loadGenericCls(readCls.ells,Cls,'ks')

        self._haveKS = True

    def loadKSDirect(self,ellsCls,Cls):
        self.theory.loadGenericCls(ellsCls,Cls,'ks')

        self._haveKS = True

    def KnoxCov(self,specTypeXY,specTypeWZ,ellBinEdges,fsky):
        '''
        returns cov(Cl_XY,Cl_WZ),signalToNoise(Cl_XY)^2, signalToNoise(Cl_WZ)^2
        '''
        def ClTot(spec,ell):
            a,b = spec
            if a==b:
                Noise = self.Nls[spec](ell)
            else:
                Noise = 0.
            return self.theory.gCl(spec,ell)+Noise        
        
        X, Y = specTypeXY
        W, Z = specTypeWZ

        ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
        ellWidths = np.diff(ellBinEdges)

        covs = []
        sigs1 = []
        sigs2 = []
        for ellMid,ellWidth in zip(ellMids,ellWidths):
            ClSum = ClTot(X+W,ellMid)*ClTot(Y+Z,ellMid)+ClTot(X+Z,ellMid)*ClTot(Y+W,ellMid)
            var = ClSum/(2.*ellMid+1.)/ellWidth/fsky
            #var = var/1.e3
            covs.append(var)
            sigs1.append(self.theory.gCl(specTypeXY,ellMid)**2./var)
            sigs2.append(self.theory.gCl(specTypeWZ,ellMid)**2./var)
        #print ellMid
        return np.array(covs), np.array(sigs1), np.array(sigs2)

    def sigmaClSquared(self,specType,ellBinEdges,fsky):
        return self.KnoxCov(specType,specType,ellBinEdges,fsky)[0]

    def sn(self,ellBinEdges,fsky,specType):
        
        var, sigs1, sigs2 = self.KnoxCov(specType,specType,ellBinEdges,fsky)

        signoise = np.sqrt(sigs1.sum())
        errs = np.sqrt(var)

        return signoise, errs
            

    def snRatio(self,ellBinEdges,fsky):
        
        ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
        ellWidths = np.diff(ellBinEdges)

        sumchisq = 0.
        signum = 0.
        sigden = 0.
        
        for ellMid,ellWidth in zip(ellMids,ellWidths):
            Clkk = self.theory.gCl('kk',ellMid)
            Nlkk = self.Nls['kk'](ellMid)
            Nlgg = self.Nls['gg'](ellMid)
            Nlss = self.Nls['ss'](ellMid)
            Clkg = self.theory.gCl('kg',ellMid)
            Clgg = self.theory.gCl('gg',ellMid)
            Clks = self.theory.gCl('ks',ellMid)
            Clss = self.theory.gCl('ss',ellMid)
            Clsg = self.theory.gCl('sg',ellMid)
    
            r0 = Clkg / Clsg
            pref = 1./(fsky*(2.*ellMid+1.)*ellWidth) # added ellWidth

            sigmaZsq = ((Clkk+Nlkk)*(Clgg+Nlgg))+(Clkg**2.)+((r0**2.)*((Clss+Nlss)*(Clgg+Nlgg)+Clsg**2.))-(2*r0*(Clks*(Clgg+Nlgg)+Clkg*Clsg))

            sigmaZsq = sigmaZsq * pref

            numer = (Clsg**2.)
            denom = sigmaZsq

            signum += (Clkg*Clsg/sigmaZsq)
            sigden += ((Clsg**2.)/sigmaZsq)


            chisq = numer/denom

            sumchisq += chisq
    
        maxlike = signum/sigden

        sigmaR = 1./np.sqrt(sumchisq)
        percentR = sigmaR*100./maxlike
        snR = maxlike/sigmaR

        return percentR,snR,maxlike
           


def getColNum(m1,m2,numwindows):
    '''
    m = 0  -> phi
    m = 1+ -> window 1, 2, etc.
    '''

    nfields = numwindows + 3

    retval = 3+nfields*(m1+2)+m2
    return retval        

def main(argv):

    def _parseSky(fstr):
        # converts e.g. 0.5f and 20000.d to fsky = 0.5

        if fstr[-1].lower()=='f':
            return float(fstr[:-1])
        elif fstr[-1].lower()=='d':
            full = 4.*np.pi*(180./np.pi)**2.            
            return float(fstr[:-1])/full
        else:
            raise ValueError
        
    def _overlapEllrange(ellranges,ellwidth,sec1,sec2):
        # determines minimal overlap ellrange
        
        lmin = max(ellranges[sec1+'-'+sec1][0],ellranges[sec2+'-'+sec2][0])
        lmax = min(ellranges[sec1+'-'+sec1][-1],ellranges[sec2+'-'+sec2][-1])

        return np.arange(lmin,lmax,ellwidth)
        

    # Load ini
    try:
        iniFile = argv[0]
    except:    
        iniFile = "input/projections.ini"
        
    config = ConfigParser()
    config.read(iniFile)
    theoryFile = config.get('theory','fileName')
    NlFile = config.get('cmb','Nlkkfile')
    bgNgal = config.getfloat('photo','ngal')
    fgNgal = config.getfloat('spec','ngal')
    shapeNoise = config.getfloat('photo','shapeNoise')
    ellwidth = config.getfloat('global','ellwidth')

    ellranges = {}
    fskys = {}

    # holds list of pairs to do S/N for
    todo = []

    # get ellranges and fskys for autos
    for sec in ['cmb','photo','spec']:
        lmin = config.getfloat(sec,'ellmin')
        lmax = config.getfloat(sec,'ellmax')
        ellranges[sec+'-'+sec] = np.arange(lmin,lmax,ellwidth)
        fskys[sec+'-'+sec] = _parseSky(config.get('theory','sky_'+sec))
        todo.append((sec,sec))

    # get ellranges and fskys for crosses
    for sec1,sec2 in zip(['cmb','cmb','photo'],['photo','spec','spec']):
        ellranges[sec1+'-'+sec2] = _overlapEllrange(ellranges,ellwidth,sec1,sec2)
        fskys[sec1+'-'+sec2] = _parseSky(config.get('theory','sky_'+sec1+'_'+sec2))
        todo.append((sec1,sec2))

    shorts = {}
    shorts['cmb'] = 'k'
    shorts['photo'] = 's'
    shorts['spec'] = 'g'
    
    LF = LensForecast()
    LF.loadKK(theoryFile,getColNum(0,0,2),NlFile,1,nnorm="none",ntranspower=[0.,1.0])
    LF.loadSS(theoryFile,getColNum(1,1,2),shapeNoise=shapeNoise,ngal=bgNgal)
    LF.loadKS(theoryFile,getColNum(0,1,2))    
    LF.loadGG(theoryFile,getColNum(2,2,2),ngal=fgNgal)
    LF.loadKG(theoryFile,getColNum(0,2,2))
    LF.loadSG(theoryFile,colnums=[getColNum(1,2,2),getColNum(2,1,2)])

    for sec1, sec2 in todo:
        short = shorts[sec1]+shorts[sec2]
        signoise = LF.sn(ellranges[sec1+'-'+sec2],fsky=fskys[sec1+'-'+sec2],specType=short)[0]
        print sec1+"-"+sec2+" S/N = ", signoise
        
        if short=='ks':
            sig = 1.
            percenterr = 100./sig**2./signoise
            print "Percent multiplicative bias = ", percenterr
        if short=='kg':
            sig = 2.
            percenterr = 100./sig**2./signoise
            print "Percent bias from CMB = ", percenterr
        if short=='sg':
            sig = 2.
            percenterr = 100./sig**2./signoise
            print "Percent bias from shear = ", percenterr
            
    

    print LF.snRatio(ellranges['photo-spec'],fskys['photo-spec'])

    
    sys.exit()
    
    # #kknoise = "/astro/astronfs01/workarea/msyriac/Planck2015/kappa_maps/nlkk.dat"
    # #kknoise = "data/s4_1uK_EB.csv"
    # kknoise = "data/advact_mv.csv"
    # ellw = 50
    # ellrange = np.arange(100,1900,ellw)+ellw/2

    # #first = "hsc"
    # first = "cfht"

    # #mid = "lowk"
    # #mid = "quick"
    # mid = "highk"

    # # fgNgal = 0.026
    # # gal = "cmass"
    # # fskyKG = 10000./42000.

    
    # fgNgal = 0.083
    # gal = "desi"
    # fskyKG = 14000./42000.

    # #fgNgal = 50.
    # #gal = "lsst"
    # #fskyKG = 20000./42000.

    # bias = 2.0

    # pl1 = Plotter()
    # pl2 = Plotter()

    # #N=6
    # #fgs = [0.9,3.7,10.9,18.8,13.2,2.2]
    # N=10
    # fgs = [5.7,15.,13.3,8.3,4.3,2.0,0.9,0.4,0.1,0.1]
    # sns = []
    # errs = []

    # for i in range(N):

    #     LF = LensForecast()

    #     fgNgal = fgs[i]
        
    #     LF.loadKK("data/planck2015_"+first+"_desi_lensing_"+mid+"_scalCovCls.dat",getColNum(0,0,2),kknoise,1,nnorm="none",ntranspower=[0.,1.0])
    #     #print "planck-planck S/N = ", LF.getSN(ellrange,fsky=0.65,specType="kk")[0]
        
    #     LF.loadGG("data/planck2015_"+first+"_"+gal+"_lensing_"+mid+"_lbin_"+str(i+1)+"_scalCovCls.dat",getColNum(2,2,2),ngal=fgNgal)
    #     LF.loadKG("data/planck2015_"+first+"_"+gal+"_lensing_"+mid+"_lbin_"+str(i+1)+"_scalCovCls.dat",getColNum(0,2,2))


    #     pl1.add(LF.ells_gg,np.array(LF.ells_gg)*np.array(LF.Clgg(LF.ells_gg)),label="bin "+str(i+1))
    #     pl2.add(LF.ells_kg,np.array(LF.ells_kg)*np.array(LF.Clkg(LF.ells_kg)),label="bin "+str(i+1))

        
    #     sn = LF.getSN(ellrange,fsky=fskyKG,specType="kg")[0]
    #     err = 1./bias/sn

    #     sns.append(sn)
    #     errs.append(err)
    #     print "planck-cmass S/N = ", sn
    #     print "error on bias = ", err

    # print sns
    # print errs
    # pl1.legendOn(loc = "upper right")
    # pl2.legendOn(loc = "upper right")
    # pl1.done("output/clggs.png")
    # pl2.done("output/clkgs.png")

    # sys.exit()

    
    ellw = 50
    ellrange = np.arange(100,1980,ellw)+ellw/2
    ellrange_act = np.arange(40,2900,ellw)+ellw/2
    ellrange_shear = np.arange(100,1500,ellw)+ellw/2

    kknl = "data/NoiseCurvesKK/act_10uk.csv"
    #kknl = "data/NoiseCurvesKK/advact_mv.csv"
    #kknl = "data/advact_mv.csv"

    
    #kknl = "data/NoiseCurvesKK/advAct_20k_white_mv.csv"
    #kknl = "data/NoiseCurvesKK/advAct_20k_realisitc_HWP_EE.csv"
    #kknl = "data/NoiseCurvesKK/advAct_20k_realisitc_NoHWP_TT.csv"

    #kknl = "data/NoiseCurvesKK/advAct_5k_HWP2x_mv.csv"
    #kknl = "data/NoiseCurvesKK/advAct_20k_realisitc_HWP_TT.csv"
    #kknl = "data/NoiseCurvesKK/advAct_20k_realisitc_NoHWP_TT.csv"
    
    fskyACT = 500./42000.
    fskyCMBShear = 500./42000.


    fskyCMB = fskyACT    
    fskySpecCMB = 500./42000.
    fskyCMASS = 500./42000.
    fskyKG = 0.2
    fgNgal = 0.026
    fgBias = 2.0
    

    # fskyShear = 1500./42000.
    # bgNgal = 9.0
    # galsurvey = "kids"
    # app = ""

        
    fskyShear = 40./42000.
    bgNgal = 15.0
    galsurvey = "hsc"
    app = ""

    # fskyShear = 20000./42000.
    # bgNgal = 30.0
    # galsurvey = "lsst"
    # app = ""

    # fsky = 140./42000.
    # bgNgal = 12.5
    # galsurvey = "cfht"

    # fskyShear = 5000./42000.
    # bgNgal = 6.0
    # galsurvey = "des"
    # app = ""


    LF = LensForecast()



    mvplanck = "/astro/astronfs01/workarea/msyriac/Planck2015/kappa_maps/nlkk.dat"

    # TT planck up to 4000
    #LF.loadKK("data/planck2015_"+galsurvey+"_cmass_lensing_scalCovCls.dat",getColNum(0,0,2),"data/Nlkk_planck_4000.csv",1,nnorm="none",ntranspower=[0.,1.0])

    # data MV planck up to 2000
    #LF.loadKK("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(0,0,2),mvplanck,1,nnorm="none",ntranspower=[0.,1.0])
    LF.loadKK("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(0,0,2),kknl,1,nnorm="none",ntranspower=[0.,1.0])
    LF.loadSS("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(1,1,2),shapeNoise=0.35,ngal=bgNgal)
    LF.loadKS("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(0,1,2))
    
    LF.loadGG("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(2,2,2),ngal=fgNgal)
    LF.loadKG("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",getColNum(0,2,2))
    LF.loadSG("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_"+app+"scalCovCls.dat",colnums=[getColNum(1,2,2),getColNum(2,1,2)])
    


    # print "planck-planck S/N = ", LF.getSN(ellrange,fsky=0.65,specType="kk")[0]
    # print "planck-shear S/N = ", LF.getSN(ellrange,fsky=fskyShear,specType="ks")[0]
    # print "planck-cmass S/N = ", LF.getSN(ellrange,fsky=fskyCMASS,specType="kg")[0]
    # print "shear-cmass S/N = ", LF.getSN(ellrange,fsky=fskyShear,specType="sg")[0]
    # print "shear-shear S/N = ", LF.getSN(ellrange_shear,fsky=fskyShear,specType="ss")[0]
    # print "cmass-cmass S/N = ", LF.getSN(ellrange,fsky=fskyCMASS,specType="gg")[0]

    print "planck-planck S/N = ", LF.sn(ellrange,fsky=fskyCMB,specType="kk")[0]
    print "planck-shear S/N = ", LF.sn(ellrange,fsky=fskyShear,specType="ks")[0]
    print "planck-cmass S/N = ", LF.sn(ellrange,fsky=fskyCMASS,specType="kg")[0]
    print "shear-cmass S/N = ", LF.sn(ellrange,fsky=fskyShear,specType="sg")[0]
    print "shear-shear S/N = ", LF.sn(ellrange_shear,fsky=fskyShear,specType="ss")[0]
    print "cmass-cmass S/N = ", LF.sn(ellrange,fsky=fskyCMASS,specType="gg")[0]

    print LF.snRatio(ellrange_shear,fskyCMBShear)

    #sys.exit()
    

    #LF.loadKK("data/planck2015_"+galsurvey+"_cmass_lensing_scalCovCls.dat",getColNum(0,0,2),"data/Nldd_act_16uK_3000.0.csv",1,nnorm="lsq",ntranspower=[2.0,0.25])

    #LF.loadKK("data/planck2015_"+galsurvey+"_cmass_lensing_scalCovCls.dat",getColNum(0,0,2),"data/act_10uk.csv",1)
    #LF.loadKK("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_scalCovCls.dat",getColNum(0,0,2),"data/NoiseCurvesKK/act_8uK_mv.csv",1)
    #LF.loadKK("data/TheorySpectra/planck2015_"+galsurvey+"_cmass_lensing_scalCovCls.dat",getColNum(0,0,2),kknl,1)



    
    print 'done'

if (__name__ == "__main__"):
    main(sys.argv[1:])
