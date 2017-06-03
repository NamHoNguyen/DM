import numpy as np
import matplotlib.pyplot as plt
import itertools

colors = itertools.cycle(['b', 'g', 'r', 'm', 'y', 'c', 'k'])

power11=np.loadtxt('../../software/HMcode/Aug1_cdm_mead1.dat',unpack=True)
#power12=np.loadtxt('../../software/HMcode/Aug1_cdm_mead2.dat',unpack=True)
#power1=np.loadtxt('output/Aug1_cdm_camb_mead_lin.dat',unpack=True)
#power21=np.loadtxt('output/Aug1_cdm_camb.dat',unpack=True)
#power22=np.loadtxt('output/Aug1_cdm_camb_mead2.dat',unpack=True)
k11=power11[0,:]
#k12=power12[0,:]
#k21=power21[0,:]
#k22=power22[0,:]
meadshape=np.shape(power11)
num_reds=meadshape[0]-1
for i in range(1,num_reds+1):
    d11=power11[i,:]
#    d12=power12[i,:]
#    d21=power21[i,:]
#    d22=power22[i,:]
    c = colors.next()
    plt.plot(k11,d11,'--'+c)
#    plt.plot(k12,d12,'--'+c)
#    plt.plot(k21,d21,'-'+c)
#    plt.plot(k22,d22,'-'+c)
plt.xscale('log')
plt.yscale('log')
plt.show()
