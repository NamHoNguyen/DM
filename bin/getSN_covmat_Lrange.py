'''
Get SNR as a function of Lmin and Lmax
Note:
- Have an option to use only the diagonal elements
'''
import numpy as np
from sn_forecast import LensingForecast

ClsCDM = np.loadtxt('data/May21_matter2lens_WF_CDM_cut_ibarrier_iconc_fCls.csv')
ClsFDM = np.loadtxt('data/May21_matter2lens_WF_FDM_1.0_cut_ibarrier_iconc_fCls.csv')
LF = LensingForecast(ClsCDM,ClsFDM)

beam = 0.3
fsky = 0.1
noise = 0.1
diag = False
saveDate = 'June9'

if diag: print '-------->Just diagonal terms!'
#print 'beam,fksy,noise,S/N'
    
ellBinEdges = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_4sqdeg_lbin_edges.npy')
ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2

cov = np.load('data/experiment_'+str(beam)+'arc_'+str(noise)+'uk_4sqdeg_covmat_rescaled.npy')
cov = cov*4./41253./fsky # scale with the right fsky

# Zero out off-diagonals
if diag:
    for i in range(len(cov)):
        for j in range(len(cov)):
            if i!=j: cov[i,j] = 0
            
# Calculate SNR squared
sn2 = LF.SNR2_fromCov(cov,ellBinEdges)

snLmin = np.zeros([len(ellMids),2])
snLmax = np.zeros([len(ellMids),2])

for i in range(len(sn2)):
    # Getting the lower right sub matrix
    sn2_sub_Lmin = sn2[i:,i:]
    snLmin[i] = [ellBinEdges[i],np.sqrt(np.sum(sn2_sub_Lmin))]

    # Getting the upper left sub matrix
    sn2_sub_Lmax = sn2[:i+1,:i+1]
    snLmax[i] = [ellBinEdges[i+1],np.sqrt(np.sum(sn2_sub_Lmax))]
    
np.savetxt('data/'+saveDate+'_SNvsLMin.csv',snLmin)
np.savetxt('data/'+saveDate+'_SNvsLMax.csv',snLmax)
#print beam,fsky,noise,ellBinEdges[i],np.sqrt(np.sum(sn2_sub))
