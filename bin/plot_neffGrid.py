import numpy as np
from driver import FisherForecast
import matplotlib.cm as cm
import matplotlib.pyplot as plt

try:
    iniFile = argv[0]
except:
    iniFile = 'input/fisher_neffGrid.ini'
prefix = ''
dClsRoot ={}
name = 'May30_neffGrid_mv_nocutS4_nonlinear'

#noiseRange = np.linspace(0.1,1.0,10)
#fskyRange = np.linspace(0.025,0.1,4)
#noiseRange = [0.1,0.3,0.5,0.8,1.0]
#fskyRange = [0.025,0.05,0.1,0.4]
noiseRange = [1.0] #,0.5]
fskyRange = [0.4] #,0.1]
#expOverride = 'DM-18arcsec'
expOverride = 'S4'

data = np.zeros([len(fskyRange),len(noiseRange)])
j = 0
for noise in noiseRange:
    i = 0
    for fsky in fskyRange:
        print '-----noise,fksy:',noise,fsky,'-----'
        nlkkLoc = {expOverride:'output/dump/May30_gradCut_2000_polComb_mv_beamY_2.0_noiseY_'+str(noise)+'_tellminY_30.0_tellmaxY_30000.0_kmax_30000.txt'}
        #nlkkLoc = {expOverride:'output/dump/May30_gradCut_2000_polComb_mv_beamY_0.3_noiseY_'+str(noise)+'_tellminY_100_tellmaxY_45000_kmax_60000.txt'}
        fsky = {expOverride:fsky,'HighEllPlanck':(0.6-fsky)}
        F = FisherForecast(iniFile,prefix=prefix,dClsRoot=dClsRoot,nlkkLocOverride=nlkkLoc,fskyOverride=fsky,noiseTOverride=noise)
        FisherMat = F.calcFisher(verbose = True)
        data[i,j]=F.margSigma('nnu')
        i+=1
    j+=1
np.savetxt("output/"+name+"_data.csv",data,delimiter=",")
print data
'''
# Color grid plot
val = np.loadtxt("output/"+name+"_data.csv",delimiter=",")
#print val
val = np.flipud(val)
#vmin,vmax=np.floor(val.min()),np.ceil(val.max())
vmin,vmax=val.min()*0.9,val.max()*1.1
print vmin,vmax
#vmin,vmax=15.0,35.0
print val
fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.imshow(val,extent=(noiseRange.min()-0.1,noiseRange.max()+0.1,fskyRange.min()-0.025,fskyRange.max()+0.025),cmap=cm.YlGnBu,interpolation='nearest',vmin=vmin,vmax=vmax,aspect='auto')
#im = ax.imshow(val,extent=None,cmap=cm.YlGnBu,interpolation='nearest',vmin=vmin,vmax=vmax)
ax.set_xlabel('Temp Noise ($\mu$K-arcmin)')
ax.set_ylabel('Fraction of sky')
ax.set_title('Constraint on $N_{\\rm eff}$')
val = np.around(np.flipud(val),3)
print val

for i in range(len(fskyRange)):
    for j in range(len(noiseRange)):
        ax.text(noiseRange[j],fskyRange[i],val[i,j],ha='center',va='center')

ax.set_xticks(noiseRange)
ax.set_yticks(fskyRange)
print fskyRange,noiseRange
fig.colorbar(im)
figname = 'output/'+name+'.png'
#fig.savefig(figname,format='png')
print "Saved fig ",figname
plt.show()                         
'''
