import numpy as np
import matplotlib.pyplot as plt
import camb
from camb import model
import os,sys
from scipy.interpolate import interp2d
import itertools

colors = itertools.cycle(['b','g','r','c','m','y','k'])
a = np.loadtxt(os.environ['AXIONCAMB_DIR']+'Apr16_CDM_matterpower_all.dat')
#z = 12. #0.06
zs = np.linspace(0.06,12.,200)
omm = 0.2993
H0 = 70.2
h = H0/100.
c = 3.e5

kmax = 30
ombh2 = 0.0225
omch2 = 0.125
ns = 0.961
tau = 0.06

ks = a[:,0]*h
Ps = a[:,1::][:,::-1]/h**3
#P_phi=9.*(omm**2)*((H0/h)**4  )/8./np.pi**2 * (1./ks)*(1.+z)**2 * Ps/c**4
P_phi=9.*(omm**2)*((H0)**4  )/4. *(1.+zs)**2/c**4 *Ps
P_phi=interp2d(zs,ks,P_phi,bounds_error=True)
print max(ks)
'''
k=np.exp(np.log(10)*np.linspace(-3.9,1.69,200))
zs = np.linspace(0.06,12.,7)
for z in zs:
    Pk1 = P_phi(z,k)
    plt.loglog(k, Pk1,label='axion')
    plt.xlim([1e-4,kmax])
    plt.xlabel('k Mpc')
    plt.ylabel('$P_\Psi\, Mpc^{-3}$')
    plt.legend(['z=%s'%z for z in zs]);
        
#plt.legend(loc='lower left')
plt.show()
sys.exit()
#plt.loglog(ks,P_phi/c**4,label='axion')
'''
pars = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2,tau=tau,TCMB=2.7252,YHe=0.24)
pars.InitPower.set_params(ns=ns,As = 2.05e-09)
PK = camb.get_matter_power_interpolator(pars, nonlinear=False,hubble_units=False, k_hunit=False, kmax=kmax,var1=model.Transfer_Weyl,var2=model.Transfer_Weyl, zmax=max(zs))


k=np.exp(np.log(10)*np.linspace(-3.9,np.log10(kmax),200))
#print k

zs = np.linspace(0.06,12.,7)
for z in zs:
    color = colors.next()
    Pk1 = P_phi(z,k)
    Pk2 = PK.P(z,k)
    #print np.mean(Pk1/Pk2),np.mean(Pk2/Pk1)
    print max(Pk1)/max(Pk2),max(Pk2)/max(Pk1)
    print min(Pk1)/min(Pk2),min(Pk2)/min(Pk1)
    plt.loglog(k, Pk1,color,label='axion')
    plt.loglog(k, Pk2,color+'--',label='pycamb')
    plt.xlim([1e-4,kmax])
    plt.xlabel('k Mpc')
    plt.ylabel('$P_\Psi\, Mpc^{-3}$')
    #plt.legend(['z=%s'%z for z in zs]);
        
plt.legend(loc='lower left',ncol=3)
plt.show()
