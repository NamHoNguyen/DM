import numpy as np
import matplotlib.pyplot as plt

jia = np.loadtxt('ell_recon_raw_autoN0subbed_jia.csv')
jia_grf = np.loadtxt('ell_recon_raw_autoN0subbed_jia_grf.csv')

plt.loglog(jia[:,0],jia[:,3],label='auto jia')
plt.loglog(jia_grf[:,0],jia_grf[:,3],label='auto jia grf')
plt.semilogx(jia_grf[:,0],abs(jia[:,3]-jia_grf[:,3]),label='abs(auto diff)')
plt.legend()
#plt.xlim([100,49000])
#plt.ylim([1.e-10,1.])
#plt.show()
plt.savefig('output/auto_diff.png')
