import numpy as np
import os,ConfigParser
from Limber import Matter2Lens_Limber
from scipy.interpolate import interp1d
def getClkkWF(params):
    
    # Make params.f90
    script = '!! This is a params file for WFCode\n!! Entries must be in the right order until I sort out a Makefile'
    script += '\n\t'+str(params['omega_m'])+'\t ! omega_m'
    script += '\n\t'+str(params['omega_b'])+'\t ! omega_b'
    script += '\n\t'+str(params['h'])+'\t ! h'
    script += '\n\t'+str(params['w'])+'\t ! w'
    script += '\n\t'+str(params['sigma8'])+'\t ! sigma8'
    script += '\n\t'+str(params['ns'])+'\t ! ns'
    script += '\n\t'+str(params['wa'])+'\t ! wa'
    script += '\n!! WDM params only used if iwdm = 1 in inputs'
    script += '\n\t'+str(params['wdm_m'])+'\t ! wdm mass in keV'
    script += '\n\t'+str(params['wdm_deg'])+'\t ! wdm degrees of freedom, =1.5 for spin 1/2 fermion'
    script += '\n!! FDM params only used if ifdm = 1 in inputs '
    script += '\n\t'+str(params['fdm_m'])+'\t ! fdm mass in 1e-22 eV'
    
    filename = os.environ['HMCODE_DIR']+'params.f90'
    with open(filename,'w') as tempFile:
        tempFile.write(script)

    filename = os.environ['HMCODE_DIR']+params['output']+'_params.f90'
    with open(filename,'w') as tempFile:
        tempFile.write(script)

    # Make inputs.f90
    script = '!! This is an inputs file for WFCode\n!! Must be in right order until I sort out a Makefile'
    script += '\n\t\''+params['output']+'.dat\'\t ! output file name '
    script += '\n\t'+str(params['nk'])+'\t ! nk number of k points, log spaced'
    script += '\n\t'+str(params['kmin'])+'\t ! kmin'
    script += '\n\t'+str(params['kmax'])+'\t ! kmax'
    script += '\n\t'+str(params['nz'])+'\t ! nz number of z points, log spaced'
    script += '\n\t'+str(params['zmin'])+'\t ! zmin'
    script += '\n\t'+str(params['zmax'])+'\t ! zmax'
    script += '\n\t'+str(params['imead'])+'\t ! imead'
    script += '\n\t'+str(params['iwdm'])+'\t ! iwdm'
    script += '\n\t'+str(params['ifdm'])+'\t ! ifdm'
    script += '\n\t'+str(params['ibarrier'])+'\t ! ibarrier'
    script += '\n\t'+str(params['iconc'])+'\t ! iconc'

    filename = os.environ['HMCODE_DIR']+'inputs.f90'
    with open(filename,'w') as tempFile:
        tempFile.write(script)

    filename = os.environ['HMCODE_DIR']+params['output']+'_inputs.f90'
    with open(filename,'w') as tempFile:
        tempFile.write(script)

    # Call WarmAndFuzzy
    os.chdir(os.environ['HMCODE_DIR'])
    os.system('gfortran -o wf HM_warm_fuzzy.f90; ./wf')
    os.chdir(os.environ['DM_DIR'])
    #os.system('gfortran -o wf '+os.environ['HMCODE_DIR']+'HM_warm_fuzzy.f90; '+os.environ['HMCODE_DIR']+'wf')

    # Call Limber
    H0 = float(params['h'])*100.
    omm = float(params['omega_m'])
    lmin = float(params['lmin'])
    lmax = float(params['lmax'])
    nl = float(params['nl'])
    log = False
    zArr = []
    whatPS = 2 # always 2 for WarmAndFuzzy Matter PS
    
    powerFile = os.environ['HMCODE_DIR']+params['output']+'.dat'
    Clkk = Matter2Lens_Limber(powerFile,H0,omm,lmin,lmax,nl,zArr,whatPS=whatPS,log=log)
    return Clkk
    
# Main Program
iniFile = 'input/makeDerivsWF.ini'
Config = ConfigParser.SafeConfigParser()
Config.optionxform = str
Config.read(iniFile)

paramList = []    
fparams = {}
stepSizes = {}

for (key, val) in Config.items('WarmAndFuzzy'):
    if ',' in val:
        param, step = val.split(',')
        paramList.append(key)
        fparams[key] = param
        stepSizes[key] = float(step)
    else:
        fparams[key] = val
        
print fparams
print paramList,stepSizes

ellRange = range(int(fparams['lmax'])+1)

# Save fiducials    
print "Calculating and saving fiducial cosmology..."
fClkk = getClkkWF(fparams)
fCls = interp1d(fClkk[:,0],fClkk[:,1],kind="linear",bounds_error=False,fill_value=0.)
np.savetxt('data/'+fparams['output']+'_WF_fCls.csv',np.vstack([ellRange,fCls(ellRange)]).T,delimiter=',')
                    
# Calculate and save derivatives
for paramName in paramList:
    h = stepSizes[paramName]
    print "Calculating forward difference for ", paramName
    pparams = fparams.copy()
    pparams[paramName] = str(float(fparams[paramName]) + 0.5*h)
    pClkk = getClkkWF(pparams)
                        
    print "Calculating backward difference for ", paramName
    mparams = fparams.copy()
    mparams[paramName] = str(float(fparams[paramName]) - 0.5*h)
    mClkk = getClkkWF(mparams)

    print pClkk,mClkk
    dClkk = (pClkk-mClkk)/h
    dClkk[:,0] = pClkk[:,0]
    dCls = interp1d(dClkk[:,0],dClkk[:,1],kind="linear",bounds_error=False,fill_value=0.)
    np.savetxt('data/'+fparams['output']+'_WF_dCls_'+paramName+'.csv',np.vstack([ellRange,dCls(ellRange)]).T,delimiter=',')
        
