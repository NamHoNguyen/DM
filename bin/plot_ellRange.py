import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 18})
mpl.rcParams['lines.linewidth'] = 3

fig = plt.figure(figsize=(8,7))

dataSum = 118.012627739 #75.921184986

data1 = np.loadtxt('data/June3_SNvsLMin.csv')
l, = plt.plot(data1[:,0],data1[:,1]/dataSum*100.,'#1d8348',label='S/N vs. $L_{\\rm min}$')

data2 = np.loadtxt('data/June3_SNvsLMax.csv')
l, = plt.plot(data2[:,0],data2[:,1]/dataSum*100.,'#4a235a',label='S/N vs. $L_{\\rm max}$')

data1 = np.loadtxt('data/June3_SNvsTMin.csv')
l, = plt.plot(data1[:,0],data1[:,1]/dataSum*100.,'#1d8348',label='S/N vs. $\ell_{\\rm min}$')
l.set_dashes([10,5])

data2 = np.loadtxt('data/June3_SNvsTMax.csv')
l, = plt.plot(data2[:,0],data2[:,1]/dataSum*100.,'#4a235a',label='S/N vs. $\ell_{\\rm max}$')
l.set_dashes([10,5])

plt.xlabel('$L$ or $\ell$',fontsize=24)
plt.ylabel('% of Signal-to-Noise',fontsize=24)
plt.legend(loc='center',bbox_to_anchor=[0.75, 0.7])
plt.xlim([0,45000])
plt.locator_params(axis='x', nbins=5)
plt.savefig('output/June3_SNvsLTRange_beam0.3_noise0.1.pdf')
#plt.show()
