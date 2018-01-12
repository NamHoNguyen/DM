import numpy as np
import camb
from camb import model, initialpower
import sys,os
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d

import itertools

colors = itertools.cycle(['b','g','r','c','m','y','k'])

def Matter2Lens_Limber(powerFile,H0,omm,lmin,lmax,nl,zArr=[],omk=0,whatPS=-1,log=True):
    c = 2.99792458e5
    '''
    # Planck 2015
    ombh2 = 0.02222 #0.0225
    omch2 = 0.1197 #0.125
    ns = 0.9655 #0.961
    tau = 0.078 #0.06
    As = 2.1955e-9
    '''
    # WMAP7 - 2nd Column - Table 1 - Komatsu et. al.
    ombh2 = 0.02253
    omch2 = 0.1122
    ns = 0.967
    tau = 0.085
    As = 2.42e-9
    
    #whatPS options:
    # 0 - potential
    # 1 - matter (axionCAMB)
    # 2 - fractional comoving density perturbation (WarmAndFuzzy)
    if whatPS == -1:
        print 'Argument whatPS is not set. 0 for potential, 1 for perturbation, 2 for matter.'
        sys.exit()
    
    # This code assumes flat universe
    #omk = 0
    #First set up parameters as usual
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2,tau=tau,TCMB=2.7252,YHe=0.24)
    pars.InitPower.set_params(ns=ns,As=As) #2.05e-09)
    results= camb.get_background(pars)
    
    # Get redshifts
    f = open(powerFile)
    zrange = f.readline().split()
    zrange = np.array([float(x) for x in zrange[1:]])
    dz = zrange[1]-zrange[0]
    f.close()
    #zrange=zrange[:49]
    
    Ps = np.loadtxt(powerFile)
    print np.shape(Ps)
    #PK = interp2d(z,ks,Ps, bounds_error=True)
    # Omega radiation
    omr = omm/(results.get_derived_params()['zeq']+1.)
    Hz = H0*np.sqrt(omr*(1.+zrange)**4 + omm*(1.+zrange)**3 + (1.-omm-omr))

    if whatPS == 1:
        print ' Input Matter PS'
        ks = Ps[:,0]*(H0/100.) # 1/Mpc
        Ps = Ps[:,1:]/(H0/100.)**3 # Mpc^3
        P_phi=9.*(omm**2)*((H0)**4  )/4. *(1.+zrange)**2/c**4 *Ps
    elif whatPS == 2:
        print ' Input Density PS'
        ks = Ps[:,0]*(H0/100.) # 1/Mpc
        Ps = Ps[:,1:]*(2.*np.pi**2)/(ks.reshape([len(ks),1]))**3 # Mpc^3
        P_phi=9.*(omm**2)*((H0)**4  )/4. *(1.+zrange)**2/c**4 *Ps
        
    P_phi=interp2d(zrange,ks,P_phi,bounds_error=True,kind='linear')
    print min(ks),max(ks),min(zrange),max(zrange)

    # Can change this number to change the L range++
    nz = 200 #number of steps to use for the radial/redshift integration
    kmax=int(max(ks))  #kmax to use


    #PK = camb.get_matter_power_interpolator(pars, nonlinear=False,hubble_units=False, k_hunit=False, kmax=kmax,var1=model.Transfer_Weyl,var2=model.Transfer_Weyl, zmax=max(zrange))
    
    #For Limber result, want integration over \chi (comoving radial distance), from 0 to chi_*.
    #so get background results to find chistar, set up arrage in chi, and calculate corresponding redshifts
    chistar = results.conformal_time(0)- model.tau_maxvis.value

    
    chis = np.linspace(0,chistar,nz)
    zs=results.redshift_at_comoving_radial_distance(chis)
    #Calculate array of delta_chi, and drop first and last points where things go singular
    dchis = (chis[2:]-chis[:-2])/2
    chis = chis[1:-1]
    zs = zs[1:-1]
    # z range selection:
    imin = np.min(np.where(zs>min(zrange)))
    imax = np.max(np.where(zs<max(zrange)))
    dchis = dchis[imin:imax]
    chis = chis[imin:imax]
    zs = zs[imin:imax]
    zs0 = zs*1.
    #Get lensing window function (flat universe)
    win = ((chistar-chis)/(chis**2*chistar))**2
    #Do integral over chi
    #ls = np.linspace(lmin,lmax, (lmax-lmin)/dl+1)
    if log:
        ls = np.logspace(np.log10(lmin),np.log10(lmax),nl)
    else:
        if lmin <= 100:
            nllow = (100-lmin+1)/2
            llow = np.logspace(np.log10(lmin),np.log10(100),nllow)
            ls = np.append(llow,np.linspace(101,lmax,nl-nllow))
        else:
            ls = np.linspace(lmin,lmax,nl)
    ellmax = max(ks*min(chis))
    #print min(chis),max(chis),min(ks),max(ks)
    print 'Allowed lmax = ',max(ks*min(chis))
    print 'Allowed lmin = ',min(ks*max(chis))
        
    #cl_kappa=np.zeros(ls.shape)
    cl_kappa=np.zeros([len(ls),len(zs)])
    w = np.ones(chis.shape) #this is just used to set to zero k values out of range of interpolation
    if len(zArr) != 0:
        index = 0
        ClsArr = np.zeros([len(ls),len(zArr)])

    test = False
    for i, l in enumerate(ls):
        sys.stdout.write("\r" + 'l = ' + str(l))
        sys.stdout.flush()
        if l > ellmax:
            chisrange = np.where(chis>(l+0.5)/max(ks))
            chis = chis[chisrange]
            dchis = dchis[chisrange]      
            zs = zs[chisrange]
            zs0 = zs*1.
            win = win[chisrange]
            w = w[chisrange]
        
        k=(l+0.5)/chis        
        w[:]=1
        w[k<1e-4]=0
        w[k>=kmax]=0
        N = len(chis)
        if test:
            sys.exit()
            #cl_kappa[i] = np.dot(dchis, w*PK.P(zs, k, grid=False)*win/k**4)
            #for x,y in zip(zs,k):
            #print '----',PK.P(x,y),P_phi(x,y)
            #print '--PK--',PK.P(zs, k,grid=False),[float(P_phi(x,y)) for x,y in zip(zs, k)]
        else:
            P_phi_data = [float(P_phi(x,y)) for x,y in zip(zs, k)]
            cl_kappa[i,(len(cl_kappa[0,:])-len(zs)):] = dchis*P_phi_data*w*win/k**4
    if test:
        cl_kappa = (ls*(ls+1))**2 * cl_kappa[:,0]
        return np.transpose(np.vstack([ls,cl_kappa]))
    '''
    print cl_kappa
    print cl_kappa[:,np.where(zs0<zmaxArr[1])],cl_kappa[:,np.where(zs0<zmaxArr[1])].sum(axis=1),cl_kappa[:,np.where(zs0<zmaxArr[1])].sum(axis=1).sum(axis=1)
    sys.exit()
    '''
    if len(zArr) != 0:
        for zmin,zmax in zArr:
            izmin = -1
            izmax = -1
            nz = len(zs0)
            #print '------------\n'
            for i in range(nz):
                if (izmin!=-1) and (izmax!=-1):
                    break
                if (zs0[i]>zmin) and (izmin==-1):
                    izmin = i
                if (zs0[nz-i-1]<zmax) and (izmax==-1):
                    izmax = nz-i
            #print zmin,zmax,zs0[izmin:izmax]
            #print ClsArr[:,index],cl_kappa[:,np.where(zs<zmax)].sum(axis=1).sum(axis=1)
            #ClsArr[:,index]=cl_kappa[:,np.where(zs0<zmax)].sum(axis=1).sum(axis=1) * (ls*(ls+1))**2
            ClsArr[:,index]=cl_kappa[:,izmin:izmax].sum(axis=1) * (ls*(ls+1))**2
            index+=1
        return np.hstack([ls.reshape(len(ls),1),ClsArr])
    else:
        cl_kappa = (ls*(ls+1))**2 * cl_kappa.sum(axis=1)
        cl_kappa = interp1d(ls,cl_kappa,bounds_error=True)
        ls=np.arange(lmin,lmax+1,1)
        return np.transpose(np.vstack([ls,cl_kappa(ls)]))

    cl_limber= 4*cl_kappa/2/np.pi #convert kappa power to [l(l+1)]^2C_phi/2pi (what cl_camb is)

def main(argv):
    
    powerFiles = ['Oct16_fdm_0.5_cut_ibarrier_iconc.dat','Oct16_fdm_1.5_cut_ibarrier_iconc.dat'] #['May21_cdm_1.0_cut_ibarrier_iconc.dat']#['May21_wdm_1.0_cut_ibarrier_iconc.dat'] #'May21_cdm_1.0_cut_ibarrier_iconc.dat','May21_fdm_1.0_cut_ibarrier_iconc.dat']  #'May12_fdm_1.0_cut_zmax14_ibarrier_iconc.dat'] #'May12_cdm_0.1_cut_zmax14.dat',
    labels  = ['WF_FDM_0.5_cut_ibarrier_iconc','WF_FDM_1.5_cut_ibarrier_iconc']#['WF_CDM_1.0_cut_ibarrier_iconc']#['WF_WDM_1.0_cut_ibarrier_iconc'] #'WF_CDM_1.0_cut_ibarrier_iconc','WF_FDM_1.0_cut_ibarrier_iconc'] #'WF_CDM_cut_zmax14',
    whatPS = [2,2,2,2]
    
    zmaxArr = [[],[],[],[],[],[],[]]
    
    # Planck 2015
    omm = 0.315 #0.300688
    H0 = 67.31 #70.2

    lmin = 10
    lmax = 60000
    '''
    powerFiles = ['powtable_AGN_WMAP7_ready.dat','powtable_DMONLY_WMAP7_ready.dat']
    labels  = ['AGN_WMAP7','DMONLY_WMAP7']
    whatPS = [1,1]

    zmaxArr = [[],[],[],[],[],[],[]]
    
    # WMAP7 - 3rd Column Table 1 - Komatsu et. al.
    omm = 0.272
    H0 = 70.4
    
    lmin = 10
    lmax = 60000 #25600
    '''
    nl = 1000
    log = False

    #zmins = [0.,2.,4.,6.,0.]
    #zmaxs = [2.,4.,6.,500.,500.]
    zmins = [0.,1.,2.,3.,4.,5.,6.]
    zmaxs = [1.,2.,3.,4.,5.,6.,500.]
    zArr = zip(zmins,zmaxs)
    zArr = []
    #zmaxArr = [[],[],[0.5,1.,2.,3.,4.,5.,10.,15.,20.,600.]] #[2.,4.,5.,10.,12.,15.,20.]
    #zmaxArr = [[],[],[0.1,5.,12.,20.,600.]] #[2.,4.,5.,10.,12.,15.,20.]
    #zmaxArr = [4.,12.]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(powerFiles)):
        #if i>0: continue
        color = colors.next()
        #if i!=2: continue
        #powerFile = os.environ['HMCODE_DIR']+powerFiles[i]
        powerFile = powerFiles[i]
        label = labels[i]
        try:
            dummyFile = np.loadtxt(os.environ['AXIONCAMB_DIR']+powerFile)
            powerFile = os.environ['AXIONCAMB_DIR']+powerFile
        except:
            dummyFile = np.loadtxt(os.environ['HMCODE_DIR']+powerFile)
            powerFile = os.environ['HMCODE_DIR']+powerFile
        Clkk = Matter2Lens_Limber(powerFile,H0,omm,lmin,lmax,nl,zArr,whatPS=whatPS[i],log=log)
        #print Clkk
        l = Clkk[:,0]
        
        if len(zArr)!=0:
            for j in range(len(zArr)):
                #ax.semilogx(l,(l*(l+1.))**2*Cls[:,j+1],label='z<'+str(zmaxArr[j]))
                #print j
                #print Clkk[:,j+1]
                ax.loglog(l,Clkk[:,j+1],label='z='+str(zArr[j]))
        else:
            #ax.loglog(l,(l*(l+1.))**2*Cls[:,1]/(2.*np.pi),label=label)
            ax.loglog(l,Clkk[:,1],color,label=label)
            #ax.plot(l,(l*(l+1.))**2*Cls[:,1]/4.,label=label)
            np.savetxt('data/Oct16_matter2lens_'+label+'_fCls.csv',Clkk)
            
    ax.set_ylabel('$[\ell(\ell+1)]^2C_\ell^{\phi\phi}/4$',size=20)
    ax.set_xlabel('$\ell$',size=20)
    plt.legend(loc='lower left')
    #plt.savefig('output/Oct10_Clkk_zrange2.png')
    plt.show()
if (__name__ == "__main__"):
    main(sys.argv[1:])
            
