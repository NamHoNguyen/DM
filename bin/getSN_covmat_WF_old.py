import numpy as np
from sn_HuDedeoVale import LensForeCastHDV
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import itertools

ClFileCDM = 'data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv'
ClFile1 = 'data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv'

ClsCDM = np.loadtxt(ClFileCDM)
Cls1 = np.loadtxt(ClFile1)

#NlFile = 'data/Aug6_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.1_tellminY_100_tellmaxY_45000_kmax_60000.txt'
#NlFile = 'data/Aug6_gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.5_tellminY_100_tellmaxY_45000_kmax_60000.txt'
#Nls = np.loadtxt(NlFile)

beams = [0.3,0.1583]
#beams = [0.3]
fskys = [0.1,0.025]
noises = [0.1,0.5]
#noises = [0.5]
diag = False
dl300 = True
rescale = False
if dl300:
    dl300_str = '_dl300'
    print '-------->Using dl300!'
else: dl300_str = ''
if rescale:
    rescale_str = '_rescaled'
    print '-------->Using rescaled!'
else: rescale_str = ''


if diag: print '-------->Just diagonal terms!'
print 'beam,fksy,noise,S/N'

for beam,fsky,noise in itertools.product(beams,fskys,noises):
    
    #ellBinEdges = np.load("data/experiment_0.3arc_0.1uk_4sqdeg_lbin_edges_old.npy")
    #ellBinEdges = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_4sqdeg_lbin_edges'+dl300_str+'.npy')
    ellBinEdges = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_2000_2.91260385523sqdeg_lbin_edges'+dl300_str+'.npy')
    #ellBinEdges = np.load('dump/sims_small_scale_clkk_g_2000_cseed_0_experiment_'+str(beam)+'arc_'+str(noise)+'uk_reconstruction_small_lbin_edges_11.6491285568sqdeg'+dl300_str+'.npy')
    #ellBinEdges = np.load('data/sims_small_scale_clkk_highres_g_2000_cseed_0_experiment_'+str(beam)+'arc_'+str(noise)+'uk_reconstruction_small_unlensed_False_lbin_edges_2.91260385523sqdeg'+dl300_str+'.npy')
    ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
    N = len(ellMids)

    #cov = np.load("data/experiment_0.3arc_0.1uk_4sqdeg_covmat_old.npy")
    #cov = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_4sqdeg_covmat'+dl300_str+rescale_str+'.npy')
    cov = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_2000_2.91260385523sqdeg_covmat'+dl300_str+rescale_str+'.npy')
    #cov = np.load('dump/sims_small_scale_clkk_g_2000_cseed_0_experiment_'+str(beam)+'arc_'+str(noise)+'uk_reconstruction_small_covmat_11.6491285568sqdeg'+dl300_str+'.npy')
    #cov = np.load('data/sims_small_scale_clkk_highres_g_2000_cseed_0_experiment_'+str(beam)+'arc_'+str(noise)+'uk_reconstruction_small_unlensed_False_covmat_2.91260385523sqdeg'+dl300_str+rescale_str+'.npy')
    #cov = np.load("data/beam_0.3_noise_0.5_area_4sqdeg_covmat.npy")

    #print cov
    # rescale with the right fsky
    #cov = cov*4./41253./fsky
    #cov = cov*11.6491285568/41253./fsky
    cov = cov*2.91260385523/41253./fsky

    # Zero out off-diagonals
    if diag:
        for i in range(len(cov)):
            for j in range(len(cov)):
                if i!=j: cov[i,j] = 0

    #print '---------------COVMAT--------------'
    #print cov

    '''
    cov = np.zeros([N,N])
    ClsCDM_func = interp1d(ClsCDM[:,0],ClsCDM[:,1],kind="linear",bounds_error=True,fill_value=0.)
    Nls_func = interp1d(Nls[:,0],Nls[:,1],kind="linear",bounds_error=True,fill_value=1.e50)
    cov[range(N),range(N)] = (Nls_func(ellMid) + ClsCDM_func(ellMid))**2
    cov =  np.nan_to_num(cov)
    '''
    
    LF = LensForeCastHDV()
    LF.loadKKDirect(Cls1[:,0],Cls1[:,1],ClsCDM[:,0],ClsCDM[:,1])
    sn2 = LF.sn2Mat(cov,ellBinEdges,fsky)
    #print sn2
    print beam,fsky,noise,np.sqrt(np.sum(sn2))
    

'''
Cls1_func = interp1d(Cls1[:,0],Cls1[:,1],kind="linear",bounds_error=True,fill_value=0.)
plt.loglog(Cls1[:,0],Cls1[:,1],'o')
plt.loglog(ellMids,Cls1_func(ellMids),'*')
plt.savefig('test.png')
'''
