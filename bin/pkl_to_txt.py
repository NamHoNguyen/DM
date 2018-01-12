import pickle
import numpy as np
data = pickle.load(open('data/sendPlotN.pkl'))
#print data[0],data[1]*2.*np.pi/4.
savedata = np.array([data[0],data[1]*2.*np.pi/4.    *2.  ]).T
#savedata = np.array([data[0],data[1]*2.*np.pi/4.]).T
#print savedata
np.savetxt('data/BlakeNoisex2.csv',savedata)
