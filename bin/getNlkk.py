from ConfigParser import SafeConfigParser 
import cPickle as pickle
import numpy as np
import argparse
from lensing_interface import lensNoise
from orphics.io import list_from_config, cprint
import orphics.cosmology as cosmo
from scipy.interpolate import interp1d

# Get the name of the experiment and lensing type from command line
parser = argparse.ArgumentParser(description='Run a Fisher test.')
parser.add_argument('expName', type=str,help='The name of the experiment in input/params.ini')
parser.add_argument('lensName',type=str,help='The name of the CMB lensing section in input/params.ini. ',default="")
#parser.add_argument('saveName', nargs='?',type=str,help='Suffix for plots ',default="")
parser.add_argument('saveName',type=str,help='Suffix for plots ',default="")

args = parser.parse_args()
expName = args.expName
lensName = args.lensName
saveName = args.saveName
print expName,lensName,saveName

TCMB = 2.7255e6

# Read config
iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)

gradCut = 2000
nIter = np.inf
#tellYRange = [(100,45000)]
tellYRange = [None]
pellYOverride = (100,45000)
bellYOverride = (2000,45000)
outDir = 'data/'+saveName
#kellmaxOverride = 60000
if nIter == 1:
    iterName = '_iterOff'
else:
    iterName = '_iterOn'

cambRoot = 'data/Aug6_highAcc_CDM'
theoryOverride = cosmo.loadTheorySpectraFromCAMB(cambRoot,unlensedEqualsLensed=False,useTotal=False,TCMB = TCMB,lpad=60000,get_dimensionless=True)

for tellYOverride in tellYRange:
    ls,Nls,ellbb,dclbb,efficiency,cc,name = lensNoise(Config,expName,lensName,theoryOverride=theoryOverride,nIter=nIter,gradCut=gradCut,tellYOverride=tellYOverride,pellYOverride=pellYOverride,bellYOverride=bellYOverride)
    kellmin,kellmax = list_from_config(Config,lensName,'Lrange')
    fnKK = cosmo.noise_pad_infinity(interp1d(ls,Nls,fill_value=np.inf,bounds_error=False),kellmin,kellmax)
    Lrange = np.arange(kellmin,kellmax,10)
    fileName = outDir + name + '.txt'
    np.savetxt(fileName,np.vstack([Lrange,fnKK(Lrange)]).T)
    cprint("Delensing efficiency: "+ str(efficiency) + " %",color="green",bold=True)
    
