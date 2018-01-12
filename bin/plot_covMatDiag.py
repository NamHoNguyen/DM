import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

NlFile = 'data/Aug6_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.5_tellminY_100_tellmaxY_45000_kmax_60000.txt'
#NlFile = 'data/May27_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_45000_kmax_60000.txt'
Nls = np.loadtxt(NlFile)
#plt.loglog(Nls[:,0],Nls[:,1],label='TT')

ClFile = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
Cls = np.loadtxt(ClFile)
l,= plt.loglog(Cls[:,0],Cls[:,1],'c--',label='CDM')

ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'
Cls1 = np.loadtxt(ClFile1)
l, = plt.loglog(Cls1[:,0],Cls1[:,1],'m',label='$10^{-22}$eV FDM')

Cls_Aug6 = np.loadtxt('data/Aug6_highAcc_CDM_lenspotentialCls.dat')
plt.loglog(Cls_Aug6[:,0],Cls_Aug6[:,5]*2.*np.pi/4.,label='Aug6 $C_L^{\kappa\kappa}$')

#Cls_Nov10 = np.loadtxt('data/Nov10_highAcc_CDM_lenspotentialCls.dat')
#plt.loglog(Cls_Nov10[:,0],Cls_Nov10[:,5]*2.*np.pi/4.,label='Nov10 $C_L^{\kappa\kappa}$')


Nls_func = interp1d(Nls[:,0],Nls[:,1],kind="linear",bounds_error=True,fill_value=0.)
Cls_func = interp1d(Cls[:,0],Cls[:,1],kind="linear",bounds_error=True,fill_value=0.)
Cls1_func = interp1d(Cls1[:,0],Cls1[:,1],kind="linear",bounds_error=True,fill_value=0.)
Cls_Aug6_func = interp1d(Cls_Aug6[:,0],Cls_Aug6[:,5]*2.*np.pi/4.,kind="linear",bounds_error=True,fill_value=0.)
#Cls_Nov10_func = interp1d(Cls_Nov10[:,0],Cls_Nov10[:,5]*2.*np.pi/4.,kind="linear",bounds_error=True,fill_value=0.)


#print Nls[:,0],Cls_Aug6[:,0]
ell = np.arange(45.,45000.,1.)
plt.loglog(ell,Nls_func(ell)+Cls_Aug6_func(ell),label='Nam $C_L^{\kappa\kappa}+N_L^{\kappa\kappa}$')
# simulation

ellBinEdges = np.load('data/experiment_0.3arc_0.5uk_4sqdeg_lbin_edges_dl300.npy')
#ellBinEdges = np.load('data/sims_small_scale_clkk_g_2000_cseed_0_experiment_0.3arc_0.1uk_reconstruction_small_lbin_edges_11.6491285568sqdeg.npy')
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
ellWidths = np.diff(ellBinEdges)

# expression ~ (1/(fsky L DeltaL)) * (Clkk+Nlkk)^2 = covmat\                
cov = np.load('data/experiment_0.3arc_0.5uk_4sqdeg_covmat_dl300.npy')
#cov = np.load('data/sims_small_scale_clkk_g_2000_cseed_0_experiment_0.3arc_0.1uk_reconstruction_small_covmat_11.6491285568sqdeg.npy')
ClkkNlkk = np.sqrt( cov.diagonal() * (4./41253.) * ellMids * ellWidths )
plt.loglog(ellMids,ClkkNlkk,label='Mat old $C_L^{\kappa\kappa}+N_L^{\kappa\kappa}$')


ellBinEdges = np.load('data/experiment_0.3arc_0.5uk_2000_2.91260385523sqdeg_lbin_edges_dl300.npy')
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
ellWidths = np.diff(ellBinEdges)

cov = np.load('data/experiment_0.3arc_0.5uk_2000_2.91260385523sqdeg_covmat_dl300.npy')
ClkkNlkk = np.sqrt( cov.diagonal() * (2.91260385523/41253.) * ellMids * ellWidths )
plt.loglog(ellMids,ClkkNlkk,label='Mat rerun $C_L^{\kappa\kappa}+N_L^{\kappa\kappa}$')



#ellBinEdges = np.load('data/sims_small_scale_clkk_g_2000_cseed_0_experiment_0.3arc_0.5uk_reconstruction_small_lbin_edges_11.6491285568sqdeg.npy')
ellBinEdges = np.load('data/sims_small_scale_clkk_highres_g_2000_cseed_0_experiment_0.3arc_0.5uk_reconstruction_small_unlensed_False_lbin_edges_2.91260385523sqdeg.npy')
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
ellWidths = np.diff(ellBinEdges)

#cov = np.load('data/sims_small_scale_clkk_g_2000_cseed_0_experiment_0.3arc_0.5uk_reconstruction_small_covmat_11.6491285568sqdeg.npy')
cov = np.load('data/sims_small_scale_clkk_highres_g_2000_cseed_0_experiment_0.3arc_0.5uk_reconstruction_small_unlensed_False_covmat_2.91260385523sqdeg.npy')
ClkkNlkk = np.sqrt( cov.diagonal() * (2.91260385523/41253.) * ellMids * ellWidths )
plt.loglog(ellMids,ClkkNlkk,label='Mat new new $C_L^{\kappa\kappa}+N_L^{\kappa\kappa}$')


#Nlkk = np.sqrt( cov.diagonal() * (4./41253.) * ellMids * ellWidths ) - Cls_func(ellMids)
#Nlkk = np.sqrt( cov.diagonal() * (4./41253.) * ellMids * ellWidths ) - Cls_Aug6_func(ellMids)
#Nlkk = np.sqrt( cov.diagonal() * (4./41253.) * ellMids * ellWidths ) - Cls_Nov10_func(ellMids)
#plt.loglog(ellMids,Nlkk,label='$N_L^{\kappa\kappa}$')


plt.legend(loc='lower left')
#plt.xlim([30,1e5])
#plt.ylim([5e-13,1e-5])
plt.xlim([2000,45000])
plt.ylim([1e-11,1e-7])
#plt.savefig('Sep28_Clkk_0.5uK.png')
plt.show()
