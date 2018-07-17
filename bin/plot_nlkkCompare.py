import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('data/June21_gradCut_2000_polComb_mv_beamY_0.167_noiseY_0.1_tellminY_100_tellmaxY_45000_pellminY_100_pellmaxY_45000_kmax_60000.txt')
b = np.loadtxt('data/June21_gradCut_2000_polComb_mv_beamY_0.167_noiseY_0.1_tellminY_100_tellmaxY_45000_pellminY_2000_pellmaxY_45000_kmax_60000.txt')
c = np.loadtxt('data/June21_gradCut_2000_polComb_mv_beamY_0.167_noiseY_0.1_tellminY_100_tellmaxY_45000_pellminY_100_pellmaxY_45000_bellminY_2000_bellmaxY_45000_kmax_60000.txt')
d = np.loadtxt('data/June21_gradCut_2000_polComb_mv_beamY_0.167_noiseY_0.1_tellminY_100_tellmaxY_45000_pellminY_100_pellmaxY_45000_bellminY_1500_bellmaxY_45000_kmax_60000.txt')

plt.loglog(a[:,0],a[:,1],label='pmin=100,bmin=100')
plt.loglog(b[:,0],b[:,1],label='pmin=2000,bmin=2000')
plt.loglog(c[:,0],c[:,1],label='pmin=100,bmin=2000')
plt.loglog(d[:,0],d[:,1],label='pmin=100,bmin=1500')

plt.legend()
#plt.xlim([100,49000])
plt.ylim([1.e-10,1.])
plt.show()
