import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
import itertools
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3
dash1 = [10,5,2,5]
#dash = [[10,5],[10,5,2,5]]
color=itertools.cycle(['k','silver'])

ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'
Cls1 = np.loadtxt(ClFile1)

ClFile = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
Cls = np.loadtxt(ClFile)

cdm = interp1d(Cls[:,0],Cls[:,1],bounds_error=True)
fdm = interp1d(Cls1[:,0],Cls1[:,1],bounds_error=True)

nbins = [15]
for nbin in nbins:

    fig = plt.figure(figsize=(9,7))
    l, = plt.plot(Cls1[:,0],Cls1[:,1],'m',label='$10^{-22}$eV FDM')
    l, = plt.plot(Cls[:,0],Cls[:,1],'c--',label='CDM')
    l.set_dashes(dash1)
    
    Lrange = np.linspace(10000,40000,nbin)
    dataFile = ['Aug6_SNvsL_noise0.1.csv','Aug6_SNvsL_noise0.5.csv']
    label = ['0.1 $\mu$K, 4000 deg$^2$','0.5 $\mu$K, 4000 deg$^2$']
    ls = ['-','--']
    shiftFactor = 0.
    
    noises = [0.1,0.5]
    
    for index in range(len(noises)):
        ellBinEdges = np.load('data/experiment_0.3arc_'+str(noises[index])+'uk_2000_2.91260385523sqdeg_lbin_edges_dl300.npy')
        ell  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
        cov = np.load('data/experiment_0.3arc_'+str(noises[index])+'uk_2000_2.91260385523sqdeg_covmat_dl300.npy')
        covDiag = cov.diagonal() * (2.91260385523/41253.) / 0.1
        
        xs = []
        ys = []
        x = 0.
        y = 0.
        j = 1
        count = 0
        for i in range(len(covDiag)):
            if Lrange[j] < ell[i]:
                if count != 0:
                    xs.append(x/count)
                    ys.append(1./np.sqrt(y))
                    count = 0
                j += 1
                x = 0.
                y = 0.
            count += 1
            x += ell[i]
            y += 1./covDiag[i]
        if count != 0:
            xs.append(x/count)
            ys.append(1./np.sqrt(y))
            
        xs = np.array(xs)
        ys = np.array(ys)
        
        #yerr = np.sqrt(cov.diagonal() * (2.91260385523/41253.) / 0.1)
        print len(ys),len(ell)
        (_,caps,eb)=plt.errorbar(xs+shiftFactor,fdm(xs+shiftFactor),yerr=ys,ls='',ecolor=color.next(),capsize=5,label=label[index])
        for cap in caps:
            cap.set_markeredgewidth(2)
        eb[0].set_linestyle(ls[index])
        shiftFactor+=300
        
        
    plt.xlabel('$L$',size=24)
    plt.ylabel('$C_L^{\kappa\kappa}$',size=24)
    plt.legend(loc='upper right',ncol=1)
    plt.xlim([10000,40000])
    plt.ylim([1.3e-11,2.2e-10])
    plt.yscale('log')
    #plt.show()
    plt.savefig('output/Sep29_errorbar_nbin'+str(nbin)+'.pdf')
    plt.close()
