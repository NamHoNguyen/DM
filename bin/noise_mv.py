import numpy as np
import matplotlib.pyplot as plt

#for noise in np.linspace(0.1,1.0,10):
for noise in [1.0,2.0,2.5]:
    print 'noise=',noise
    #prefix = 'data/Aug6_gradCut_2000_polComb_'
    prefix = 'dump/Oct20_gradCut_2000_polComb_'
    suffix = '_beamY_2.0_noiseY_'+str(noise)+'_tellminY_30.0_tellmaxY_4000.0_kmax_5000.txt'
    #suffix = '_beamY_0.3_noiseY_'+str(noise)+'_tellminY_100_tellmaxY_45000_kmax_60000.txt'
    
    data = np.loadtxt(prefix+'TT'+suffix)
    data[:,1] = 0.*data[:,1]
    for est in ['TT','EE','ET','TB','EB']:
        Nls=np.loadtxt(prefix+est+suffix)
        print Nls[:,0]
        data[:,1]+= 1./Nls[:,1]
        plt.loglog(Nls[:,0],Nls[:,1],label=est)
    data[:,1]=1./data[:,1]
    for i in range(len(data[:,1])):
        if data[i,1]!=data[i,1]:
            data[i,1] = 1.e40
    plt.loglog(data[:,0],data[:,1],label='mv')
    plt.legend()
    plt.show()
    np.savetxt(prefix+'mv'+suffix,data)
