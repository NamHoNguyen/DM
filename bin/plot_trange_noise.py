# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['lines.linewidth'] = 3
import itertools

#color = itertools.cycle(['b','g','r','k'])
#ls = itertools.cycle([':','-.','--','-'])
color = itertools.cycle(['k','b','g','r'])
ls = itertools.cycle(['-','--','-.',':'])
#for x,y in [(x,y) for x in np.hstack([[1000],np.arange(5000,26000,10000)]) for y in [30000]]:
for x,y in [(x,y) for x in [100] for y in [5000,10000,20000,30000]]:
    a=np.loadtxt('output/dump/Feb19__gradCut_2000_polComb_TT_beamY_0.3_noiseY_0.1_tellminY_'+str(x)+'_tellmaxY_'+str(y)+'_kmax_30000.txt')
    for i in range(len(a[:,1])):
        if a[i,1] > 1.e40: a[i,1] =1.e40
    plt.loglog(a[:,0],a[:,1],color.next()+ls.next(),label='$\ell$ = ['+str(x)+','+str(y)+']')

plt.ylim([1.e-11,1.e-5])
plt.xlabel('$L$',size=20)
plt.ylabel('$N_L^{\kappa\kappa}$',size=20)
plt.legend(loc='lower left')
#plt.savefig('output/Mar14_tmintest.pdf')
#plt.savefig('output/Mar14_tmaxtest.pdf')
plt.savefig('output/Mar14_tmaxtest.png',dpi=200)
#plt.show()
