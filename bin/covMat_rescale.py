'''
Script rescale the diagonal elements to Nlkk that we trust. Assuming the correlation coefficients are correct
'''
import numpy as np
import itertools
from scipy.interpolate import interp1d

# interpolation is not the best way, do Mat's binning
ClsCDM = np.loadtxt('data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv')
ClsCDM_func = interp1d(ClsCDM[:,0],ClsCDM[:,1],kind="linear",bounds_error=True,fill_value=0.)


beams = [0.3,0.1583]
noises = [0.1,0.5]
fsky = 4./41253.
dl300 = True

if dl300: dl300_str = '_dl300'
else: dl300_str = ''

for beam,noise in itertools.product(beams,noises):
    print beam,noise
    name = 'data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_4sqdeg'
    ellBinEdges = np.load(name+'_lbin_edges'+dl300_str+'.npy')
    cov_old = np.load(name+'_covmat'+dl300_str+'.npy')
    Nls = np.loadtxt('data/Aug6_gradCut_2000_polComb_TT_beamY_'+str(round(beam,2))+'_noiseY_'+str(noise)+'_tellminY_100_tellmaxY_45000_kmax_60000.txt')
    Nls_func = interp1d(Nls[:,0],Nls[:,1],kind="linear",bounds_error=True,fill_value=1.e50)

    print cov_old[0,0]
    
    ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
    ellWidths = np.diff(ellBinEdges)
    
    # this is stddev/(fsky L dL)
    stddev_old = np.sqrt(cov_old.diagonal())
    stddev_new = (ClsCDM_func(ellMids) + Nls_func(ellMids)) / np.sqrt(fsky*ellMids*ellWidths)

    #print ClsCDM_func(ellMids)/Nls_func(ellMids)
    ratio = stddev_new/stddev_old
    cov_new = ratio.reshape([len(ratio),1]) * cov_old * ratio
    print cov_new[0,0]
    np.save(name+'_covmat'+dl300_str+'_rescaled.npy',cov_new)
