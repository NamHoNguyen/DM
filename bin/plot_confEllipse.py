import os
import itertools
import numpy as np
import sys
import ConfigParser
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
#import matplotlib
import matplotlib as mpl
mpl.rcParams.update({'font.family':'serif'})
mpl.rcParams.update({'font.size': 16})
mpl.rcParams['lines.linewidth'] = 3

#colors = itertools.cycle(['b', 'r', 'g', 'm', 'y', 'c', 'k'])
colors = itertools.cycle(['k', 'k', 'k', 'k', 'r'])
styles = itertools.cycle(['-','--','-.',':'])
#dashes = itertools.cycle([[1,0],[10,5,2,5]])

#matplotlib.rcParams['mathtext.default'] = 'regular'
labels = {'H0':'$H_0$','ombh2':'$\Omega_b h^2$','omch2':'$\Omega_c h^2$','ns':'$n_s$','As':'$A_s$','tau':'$\\tau$','mnu':'$\Sigma m_{\\nu}$ (eV)','nnu':'$N_{eff}$','r':'$r$'}
CL = {1:'68%',2:'95%',3:'99%'}
fontsize = 18

dataFile1 = os.environ['FISHER_DIR']+'/output/Mar31_S4_0.4_1.0_CMB+BAO_confEllipse_mnu_nnu_1sigma.csv'
dataFile2 = os.environ['FISHER_DIR']+'/output/Mar31_CSST_0.025_0.5_CMB+BAO_confEllipse_mnu_nnu_1sigma.csv'
dataFile3 = os.environ['FISHER_DIR']+'/output/Mar31_CSST_0.025_0.1_CMB+BAO_confEllipse_mnu_nnu_1sigma.csv'
dataFile4 = os.environ['FISHER_DIR']+'/output/Mar31_CSST_0.1_0.5_CMB+BAO_confEllipse_mnu_nnu_1sigma.csv'
dataFile5 = os.environ['FISHER_DIR']+'/output/Mar31_CSST_0.1_0.1_CMB+BAO_confEllipse_mnu_nnu_1sigma.csv'

dataFiles = [dataFile1,dataFile2,dataFile3,dataFile4,dataFile5]

#dataLabels = ['S4','CSST(1000,0.5)','CSST(1000,0.1)','CSST(4000,0.5)','CSST(4000,0.1)'] #,'axion unfix sep24','axion unfix sep28','axion fix sep28']
dataLabels = ['S4','1000 deg$^2$, 0.5)','(1000,0.1)','(4000,0.5)','(4000,0.1)']
alpha = {1:1.52, 2:2.48, 3:3.41}

#fig = plt.figure(figsize=(10,6))
fig = plt.figure()
ax = fig.add_subplot(111)
#param1 = 'ns'
#param2 = 'r'
 
for i in range(len(dataFiles)):
    dataFile = dataFiles[i]
    dataLabel = dataLabels[i]
    print dataFile

    config = ConfigParser.SafeConfigParser()
    config.optionxform = str
    config.read(dataFile)
  
    sigmax2 = config.getfloat('confEllipse','sigmax2')
    sigmay2 = config.getfloat('confEllipse','sigmay2')
    sigmaxy = config.getfloat('confEllipse','sigmaxy')
    try:
        param1,param2 = config.get('confEllipse','params').split(',')
    except:
        param1 = 'mnu'
        param2 = 'nnu'
    confLevel = config.getint('confEllipse','confLevel')
    angle = config.getfloat('confEllipse','angle')
    xcenter = config.getfloat('confEllipse','xcenter')
    ycenter = config.getfloat('confEllipse','ycenter')
    width = config.getfloat('confEllipse','width')
    height = config.getfloat('confEllipse','height')
    print xcenter,ycenter
    #xcenter = 0.09921434062740242
    #ycenter = 3.0389654547386447
    #xcenter = 0.06
    #ycenter = 0.12
    #ycenter = 0.08
    e = Ellipse((xcenter,ycenter),width,height, angle=angle,fill=False,label=dataLabel,linewidth=3,ls=styles.next(),color=colors.next())
    #e = Ellipse((xcenter,ycenter),width,height, angle=angle,fill=False,label=dataLabel,linewidth=3,ls='--',color=colors.next())
    
    ax.add_patch(e)
    #e.set_dashes(dashes.next())
    #e.set_hatch('/')

    ax.plot(xcenter,ycenter,'r*')#,markersize=16)

    ax.set_xlabel(labels[param1],fontsize=20)
    ax.set_ylabel(labels[param2],fontsize=20)

#ax.set_yscale('log')
#ax.set_ylim([0,4e-3])
#ax.set_xlim([0.96,0.975])
ax.set_xlim([0,0.2])
ax.set_ylim([2.9,3.2])
plt.grid()
plt.legend(loc='lower right')
#ax.set_title('Joint constraint ('+CL[confLevel]+' CL) on '+labels[param1]+' and '+labels[param2],fontsize=fontsize)
fileName = 'output/Mar31_confEllipse_'+param1+'_'+param2+'_'+str(confLevel)+'sigma'
#fileName = os.environ['FISHER_DIR']+'/output/June29_Das_confEllipse_omL_w_1sigma'
#plt.show()
#plt.savefig('Sep28_fixKT',format='png')
#plt.savefig(fileName+'.pdf')
plt.savefig(fileName+'.png',dpi=200)

#print('Saved file '+fileName+'.png')
            
