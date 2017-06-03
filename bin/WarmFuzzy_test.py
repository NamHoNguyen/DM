import sys,os
import camb
from camb import model
import numpy as np
import matplotlib.pyplot as plt
import itertools
from scipy.interpolate import interp1d

def getX(zmax,H0,omm,omk,dXdz=False,numpts=1.e5):
    c = 3.e5 # km/s
    #print zmax
    zrange = np.linspace(0.,zmax,numpts)
    #print zrange
    dz = zrange[1]-zrange[0]
    X = 0.
    for z in zrange:
        X += dz/np.sqrt( omm*(1.+z)**3 + omk*(1.+z)**2 + (1-omm-omk) )
    X = X*c/H0 # Mpc
    if dXdz:
        return X,c/H0/np.sqrt( omm*(1.+zmax)**3 + omk*(1.+zmax)**2 + (1-omm-omk) )
    else:
        return X
def getWFCls(powerFile,H0,omm,lmin,lmax,dl,zmaxArr=[]):
    # Input P_delta for linspace of z
    # This code assumes flat universe
    omk = 0
    
    f = open(os.environ['HMCODE_DIR']+powerFile)
    z = f.readline().split()
    z = [float(x) for x in z[1:]]
    dz = z[1]-z[0]
    f.close()
    
    P_delta = np.loadtxt(os.environ['HMCODE_DIR']+powerFile)
    ks = P_delta[:,0] # h/Mpc
    P_delta = P_delta[:,1:] # dimensionless
    P_phi = np.zeros(P_delta.shape)
    # Convert from density PS to potential PS
    for i in range(len(ks)):
        for j in range(len(z)):
            #P_phi[i,j] = 9.*(omm**2)*(H0**4)*((1.+z[j])**2)/4./(ks[i]**4)*P_delta[i,j]
            P_phi[i,j] = 9.*(omm**2)*(100.**4)*((1.+z[j])**2)/4./(ks[i]**4)*P_delta[i,j]
    if z[0] == 0.:
        print 'Eliminate z=0 data (first column)'
        z = z[1:]
        P_phi = P_phi[:,1:]

    # Correct the unit
    P_phi = P_phi/(3.e5)**4

    XCMB = getX(1100,H0,omm,omk)
    ells = np.linspace(lmin,lmax,(lmax-lmin)/dl+1)
    Cls = np.zeros(len(ells))
    if len(zmaxArr) != 0:
        index = 0
        ClsArr = ells.reshape(len(ells),1)

    # Print out allowed lmin and lmax
    for j in range(len(z)):
        if j == 0:
            X = getX(z[j],H0,omm,omk)
            print 'Allowed lmax =',max(ks*X)
        elif j == (len(z)-1):
            X = getX(z[j],H0,omm,omk)
            print 'Allowed lmin =',min(ks*X)
        else:
            continue
        
    for j in range(len(z)):
        if (len(zmaxArr)!=0) and (z[j]>=zmaxArr[index]):
            ClsArr = np.hstack([ClsArr,Cls.reshape(len(ells),1)])
            index += 1
        print 'Calculating for z =',z[j]
        X,dXdz = getX(z[j],H0,omm,omk,dXdz=True)
        Pz = interp1d(ks,P_phi[:,j],bounds_error=True)
        for i in range(len(ells)):
            #print 'Calculating for ell =',ells[i]
            k = ells[i]/X
            try:
                Cls[i] += (8.*(np.pi**2)/ells[i]**3)*dz*Pz(k)*X*dXdz*( (XCMB-X)/(XCMB*X) )**2
            except ValueError:
                #print 'Allowed krange: ',ks
                #print 'Input k and ell: ',k,ells[i]
                print 'Value Error. Change lmin and lmax!'
                sys.exit()
    '''
    try:
        if (len(zmaxArr)!=0) and (z[j]>=zmaxArr[index]):
            ClsArr = np.hstack([ClsArr,Cls.reshape(len(ells),1)])
    except:
        print 'Error'
    '''
    #return np.hstack([k.reshape(len(k),1),P_phi])
    #return np.hstack([k.reshape(len(k),1),P_delta])
    if len(zmaxArr)!=0:
        return ClsArr
    else:
        return np.hstack([ells.reshape(len(ells),1),Cls.reshape(len(ells),1)])

def main(argv):
    colors = itertools.cycle(['b', 'g', 'r', 'm', 'y', 'c', 'k'])
    #powerFiles = ['Aug8_cdm.dat','Aug8_fdm.dat','Aug9_cdm.dat','Aug9_fdm.dat']
    powerFiles = ['Aug8_cdm.dat','Aug8_fdm.dat','Aug9_cdm_sig8.dat','Aug9_fdm_sig8.dat']
    labels = ['CDM','FDM 1e-24eV','FDM noMod1','FDM noMod2']
    omm = 0.298 #0.275
    H0 = 69. #70.2
    lmin = 11 #30
    #lmax = 3000
    lmax = 1.e4
    dl = 1
    zmaxArr = [] #[0.5,1.,2.,3.,4.,5.,10.,12.]

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
    '''
    for name in ['Renee_cdm','Renee_fdm']:
        Renee = np.loadtxt('output/'+name+'.csv',delimiter=',')
        ax.plot(Renee[:,0],Renee[:,1]*1.e-7,label=name)
 
    for i in range(len(powerFiles)):
        powerFile = powerFiles[i]
        label = labels[i]
        Cls = getWFCls(powerFile,H0,omm,lmin,lmax,dl,zmaxArr)
        l = Cls[:,0]
        if len(zmaxArr)!=0:
            for j in range(len(zmaxArr)):
                ax.plot(l,(l*(l+1.))**2*Cls[:,j+1],label='z<'+str(zmaxArr[j]))
        else:
            ax.plot(l,(l*(l+1.))**2*Cls[:,1]/(2.*np.pi),label=label)
    '''
    Renee = np.loadtxt('output/Renee_diff.csv',delimiter=',')
    ax.plot(Renee[:,0],Renee[:,1],label="Renee's")

    axion1 = np.loadtxt('/home/nhnguyen/MEGA/software/axionCAMB/Aug10_ver11_cdm_lenspotentialCls.dat')
    axion2 = np.loadtxt('/home/nhnguyen/MEGA/software/axionCAMB/Aug10_ver11_fdm_lenspotentialCls.dat')
    axion3 = np.loadtxt('/home/nhnguyen/MEGA/software/axionCAMB/Aug10_ver12_cdm_lenspotentialCls.dat')
    axion4 = np.loadtxt('/home/nhnguyen/MEGA/software/axionCAMB/Aug10_ver12_fdm_lenspotentialCls.dat')
    ax.plot(axion1[:,0],(axion2[:,5]-axion1[:,5])/axion1[:,5],label='axionCAMB ver11')
    ax.plot(axion3[:,0],(axion4[:,5]-axion3[:,5])/axion3[:,5],label='axionCAMB ver12')
    '''
    #Cls1 = getWFCls(powerFiles[0],H0,omm,lmin,lmax,dl,zmaxArr)
    #Cls2 = getWFCls(powerFiles[1],H0,omm,lmin,lmax,dl,zmaxArr)
    Cls3 = getWFCls(powerFiles[2],H0,omm,lmin,lmax,dl,zmaxArr)
    Cls4 = getWFCls(powerFiles[3],H0,omm,lmin,lmax,dl,zmaxArr)
    l = Cls3[:,0]
    #ax.plot(l,(Cls2[:,1]-Cls1[:,1])/Cls1[:,1],label='WF zmax=12')
    ax.plot(l,(Cls4[:,1]-Cls3[:,1])/Cls3[:,1],label='WF zmax=15')
    '''
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend(loc='lower left')
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
