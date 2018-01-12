import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
from Limber import Matter2Lens_Limber
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
#mpl.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

fig = plt.figure(figsize=(10,8))
dash = [10,5]

powerFile = 'May21_cdm_1.0_cut_ibarrier_iconc.dat'
label  = 'WF_CDM_1.0_cut_ibarrier_iconc'
whatPS = 2
# Planck 2015
omm = 0.315
H0 = 67.31

lmin = 10
lmax = 60000

nl = 1000
log = False
zmins = [0,1,2,3,4,5,0]
zmaxs = [1,2,3,4,5,500,500]
zArr = zip(zmins,zmaxs)


#colors = itertools.cycle(['b','g','r','c','m','y','k'])
#color = colors.next()

powerFile = os.environ['HMCODE_DIR']+powerFile
Clkk = Matter2Lens_Limber(powerFile,H0,omm,lmin,lmax,nl,zArr,whatPS=whatPS,log=log)
l = Clkk[:,0]

for j in range(len(zArr)):
    plt.loglog(l,Clkk[:,j+1],label='z='+str(zArr[j]))
    #plt.loglog(l,Clkk[:,j+1],label='$\mathbf{z\in}$'+str(zArr[j]))
    #plt.loglog(l,Clkk[:,j+1],label='z='+str(zmins[j]))+'-'+str(zmaxs[j])

#l, = plt.plot(x[lmin:],y1,'navy',label='Linear FDM')
#l.set_dashes(dash)
                
plt.xlabel('$L$',fontsize=24)
#plt.ylabel('$[L(L+1)]^2C_L^{\phi\phi}/4$',fontsize=24)
plt.ylabel('$C_L^{\kappa\kappa,\\rm CDM}$',fontsize=24)
plt.xlim([1e2,1e5])
plt.ylim([1e-13,1e-5])
plt.legend(loc='upper right',ncol=1)
#plt.show()
plt.savefig('output/Sep25_power_zrange.pdf')    
