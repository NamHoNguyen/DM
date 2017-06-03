import sys,os
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM as FLCDM
from astropy import units as u
from astropy import constants as const
from scipy.special import spherical_jn as sphr_j

def getX(zmax,H0,omm,omk,dXdz=False,numpts=1.e5):
    c = 2.99792458e5 # km/s
    #print zmax
    '''
    zrange = np.linspace(0.,zmax,numpts)
    dz = zrange[1]-zrange[0]
    #print zrange
    X = 0.
    for z in zrange:
        X += dz/np.sqrt( omm*(1.+z)**3 + omk*(1.+z)**2 + (1-omm-omk) )
    '''
    try:
        z = np.linspace(0.,zmax,numpts)
        #dz = z[1]-z[0]
    except:
        z = np.outer(np.linspace(0.,1.,numpts),zmax)
    dz = z[1]-z[0]
    X = np.sum(dz/np.sqrt( omm*(1.+z)**3 + omk*(1.+z)**2 + (1-omm-omk) ),axis=0)
    #print X,X2,X-X2
    
    #X = X*c/H0 # Mpc
    X = (X*c/H0)*(H0/100.)  # Mpc/h
    if dXdz:
        #return X,c/H0/np.sqrt( omm*(1.+zmax)**3 + omk*(1.+zmax)**2 + (1-omm-omk) )
        return X,c/H0/np.sqrt( omm*(1.+zmax)**3 + omk*(1.+zmax)**2 + (1-omm-omk) )*(H0/100.)
    else:
        return X
    
def PS2Lens(powerFile,H0,omm,lmin,lmax,dl,zmaxArr=[],omk=0,whatPS=-1):
    c = const.c.to(u.km/u.s) # km/s
    H0 = H0 * (u.km/u.s/u.Mpc)
    #H = 100.* (u.h*u.km/u.s/u.Mpc)
    h = H0.value/100./u.h
    print h
    #whatPS options:
    # 0 - potential
    # 1 - fractional comoving density perturbation
    # 2 - matter
    if whatPS == -1:
        print 'Argument whatPS is not set. 0 for potential, 1 for perturbation, 2 for matter.'
        sys.exit()

    # This code assumes flat universe
    #omk = 0

    # Get redshifts
    f = open(powerFile)
    z = f.readline().split()
    z = np.array([float(x) for x in z[1:]])
    dz = z[1]-z[0]
    f.close()
    
    Ps = np.loadtxt(powerFile)
    ks = Ps[:,0] * (u.h/u.Mpc)# h/Mpc
    Ps = Ps[:,1:] #
    #print len(Ps),len(z)
    #print Ps,z,Ps*z
    if omk == 0:
        print "Assuming omk=0, i.e. flat universe"
        if whatPS == 0:
            print 'Input potential PS'
            P_phi = Ps * (u.km/u.s)**4
        elif whatPS == 1:
            print 'Input perturbation PS'
            P_phi = 9.*(omm**2)*(100.**4)/4.          * np.outer((1./ks**4),(1.+z)**2)/(u.h/u.Mpc)**4 * Ps
        elif whatPS == 2:
            print 'Input matter PS'
            #P_phi = 9.*(omm**2)*(100**4  )/8./np.pi**2 * np.outer( (1./ks)    , (1.+np.reshape(z,[len(z),1]))**2 ) * Ps
            P_phi = 9.*(omm**2)*((H0/h)**4  )/8./np.pi**2 * np.outer((1./ks)   ,(1.+z)**2)/(u.h/u.Mpc) * Ps*(u.Mpc/u.h)**3
            #P_phi = 9.*(omm**2)*((H0)**4  )/8./np.pi**2 * np.outer(ks**3   ,(1.+z)**2)*(u.h/u.Mpc)**3 * Ps*(u.Mpc/u.h)**3
    else:
        print "Non-flat universe is not implemented"    
    if z[0] == 0.:
        print 'Eliminate z=0 data (first column)'
        z = z[1:]
        P_phi = P_phi[:,1:]

    # Correct the unit
    P_phi = P_phi/c**4
    cosmo = FLCDM(H0=H0,Om0=omm)
    XCMB = cosmo.comoving_distance(1100)*h
    X = cosmo.comoving_distance(z)*h
    print cosmo.comoving_distance(2)
    sys.exit()
    dXdz = cosmo.hubble_distance/cosmo.efunc(z)*h

    print XCMB
    ells = np.linspace(lmin,lmax,(lmax-lmin)/dl+1)
    Cls = np.zeros(len(ells))
    if len(zmaxArr) != 0:
        index = 0
        ClsArr = ells.reshape(len(ells),1)

    # Print out allowed lmin and lmax
    #print 'Allowed lmax = ',max(ks*getX(min(z),H0,omm,omk))
    #print 'Allowed lmin = ',min(ks*getX(max(z),H0,omm,omk))
    print 'Allowed lmax = ',max(ks*min(X))
    print 'Allowed lmin = ',min(ks*max(X))

    print 'Calculating for z =',z
            
    for j in range(len(z)):
        #if (len(zmaxArr)!=0) and (z[j]>=zmaxArr[index]):
        if (len(zmaxArr)!=0) and (z[j]>=zmaxArr[index+1]):
            ClsArr = np.hstack([ClsArr,Cls.reshape(len(ells),1)])
            index += 1
        #print 'Calculating for z =',z[j]
        #X,dXdz = getX(z[j],H0,omm,omk,dXdz=True)
        Pz = interp1d(ks,P_phi[:,j],bounds_error=True)
        k = ells/X[j]
        #print ks,k
        Cls += (8.*(np.pi**2)/ells**3)*dz*Pz(k)*X[j]*dXdz[j]*( (XCMB-X[j])/(XCMB*X[j]) )**2
        '''
        for i in range(len(ells)):
            #print 'Calculating for ell =',ells[i]
            k = ells[i]/X
            #print ells[i],X[j],k
            try:
                Cls[i] += (8.*(np.pi**2)/ells[i]**3)*dz*Pz(k)*X*dXdz*( (XCMB-X)/(XCMB*X) )**2
            except ValueError:
                #print 'Allowed krange: ',ks
                #print 'Input k and ell: ',k,ells[i]
                print 'Value Error. Change lmin and lmax!'
                sys.exit()
        '''
    print Cls
    
    try:
        if (len(zmaxArr)!=0) and (z[j]>=zmaxArr[index]):
            ClsArr = np.hstack([ClsArr,Cls.reshape(len(ells),1)])
    except:
        print 'Error'
    
    #return np.hstack([k.reshape(len(k),1),P_phi])
    #return np.hstack([k.reshape(len(k),1),P_delta])
    if len(zmaxArr)!=0:
        return ClsArr
    else:
        return np.hstack([ells.reshape(len(ells),1),Cls.reshape(len(ells),1)])

def main(argv):
    '''
    H0 = 70
    omm = 0.3
    omk = 0
    z = 1
    print getX(z,H0,omm,omk)
    sys.exit()
    '''
    colors = itertools.cycle(['b', 'g', 'r', 'm', 'y', 'c', 'k'])
    #powerFiles = ['Aug8_cdm.dat','Aug8_fdm.dat','Aug9_cdm.dat','Aug9_fdm.dat']
    #powerFiles = ['Aug8_cdm.dat']
    powerFiles = [os.environ['AXIONCAMB_DIR']+'Apr16_CDM_zmax20_matterpower_all.dat']
    #powerFiles = [os.environ['AXIONCAMB_DIR']+'Apr16_CDM_zmax12_matterpower_all.dat']
    whatPS = 2
    labels = ['CDM']
    #omm = 0.298
    omm = 0.300688
    #H0 = 69.
    H0 = 70.2
    #lmin = 11 #30
    #lmax = 3000
    #lmax = 1.e4
    lmin = 2
    #lmax = 45e3
    lmax = 18000
    dl = 1
    #zmaxArr = [0.5,1.,2.,3.,4.,5.,10.,12.]
    zmaxArr = [] #[2.,4.,5.,10.,12.,15.,20.]
    #zmaxArr = [4.,12.]

    # power = getWFCls(powerFile,omm,H0)
    # k = power[:,0]
    # for i in range(1,len(power[0,:])):
    #     c = colors.next()
    #     P = interp1d(k,power[:,i],bounds_error=False,fill_value=1.e40)
    #     kk = np.logspace(np.log10(min(k)),np.log10(max(k)*0.999),5000)
    #     print k,kk
    #     plt.loglog(kk,P(kk),c+'--')
    #     plt.loglog(k,power[:,i],c+'-')
    # sys.exit()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    a=np.loadtxt(os.environ['AXIONCAMB_DIR']+'Apr10_CDM_lenspotentialCls.dat')
    #ax.plot(a[:,0],a[:,5]*2.*np.pi/4.,label='axionCls')
    ax.loglog(a[:,0],a[:,5]*2.*np.pi/4.,label='axionCls')
    #ax.semilogx(a[:,0],a[:,5]*2.*np.pi,label='axion $l^4C_l^{\phi\phi}$')
     
    for i in range(len(powerFiles)):
        #powerFile = os.environ['HMCODE_DIR']+powerFiles[i]
        powerFile = powerFiles[i]
        label = labels[i]
        Cls = PS2Lens(powerFile,H0,omm,lmin,lmax,dl,zmaxArr,whatPS=whatPS)
        l = Cls[:,0]
        if len(zmaxArr)!=0:
            for j in range(len(zmaxArr)):
                ax.loglog(l,(l*(l+1.))**2*Cls[:,j+1],label='z<'+str(zmaxArr[j]))
                #ax.semilogx(l,(l*(l+1.))**2*Cls[:,j+1],label='z<'+str(zmaxArr[j]))
        else:
            #ax.loglog(l,(l*(l+1.))**2*Cls[:,1]/(2.*np.pi),label=label)
            ax.loglog(l,(l*(l+1.))**2*Cls[:,1],label=label)
            #ax.plot(l,(l*(l+1.))**2*Cls[:,1]/4.,label=label)

    plt.legend(loc='upper right')
    plt.show()
    
if (__name__ == "__main__"):
    main(sys.argv[1:])
            
'''
# Params
H0 = 70.2
h = H0/100.
ombh2 = 0.0458*(h**2)
omch2 = (0.275-0.0458)*(h**2)
w = -1.
ns = 0.968

kmax = 1.e4#/h
kmin = 1.e-1
z = np.linspace(0.,12.,7)

pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
pars.set_dark_energy(w=w) #re-set defaults
pars.InitPower.set_params(ns=ns)
#Not non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=z, kmax=kmax)

camb.set_halofit_version(version='mead')

#Linear spectra
# pars.NonLinear = model.NonLinear_none
# results = camb.get_results(pars)
# kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
# s8 = np.array(results.get_sigma8())

#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
#pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
#results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = 2000)

N = len(kh_nonlin)
data = kh_nonlin.reshape(N,1)

print results.get_sigma8(),z_nonlin
for i in range(len(z_nonlin)):
    del2m = 4.*np.pi*(kh_nonlin/2./np.pi)**3*pk_nonlin[i,:]
    data = np.hstack([data,del2m.reshape(N,1)])
print data
np.savetxt('output/Aug1_cdm_camb_mead2.dat',data)
'''
