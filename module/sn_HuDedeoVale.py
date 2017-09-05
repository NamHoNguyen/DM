from lensingVar import TheorySpectra
import numpy as np

class LensForeCastHDV:
    def __init__(self):
        '''
        Calculate S/N for CMB lensing with non-diagonal Covariance Matrix
        '''
        self._haveKK = False   
        self.theoryTest = TheorySpectra()
        self.theoryTrue = TheorySpectra()
        self.Nls = {}

    def loadKKDirect(self,ellsClsTest,ClsTest,ellsClsTrue,ClsTrue): #,ellsNls,Nls):
        #self.Nls['kk'] = interp1d(ellsNls,Nls,bounds_error=False,fill_value=np.inf)
        self.theoryTest.loadGenericCls(ellsClsTest,ClsTest,'kk')
        self.theoryTrue.loadGenericCls(ellsClsTrue,ClsTrue,'kk')
        self._haveKK = True
                
    def sn2Mat(self,cov,ellBinEdges,fsky):
        '''
        (S/N)^2 matrix in L and L'
        '''
        N = len(cov[:,0])
        invcov = np.linalg.inv(cov)
        ellMids  =  (ellBinEdges[1:] + ellBinEdges[:-1]) / 2
        ellWidths = np.diff(ellBinEdges)
        #sn2 = 0.
        spec = 'kk'
        
        for i in range(N):
            for j in range(N):
                invcov(i,j) = (self.theoryTest.gCl(spec,ellMids[i]) - self.theoryTrue.gCl(spec,ellMids[i])) * ellWidths(i)\
                              * (self.theoryTest.gCl(spec,ellMids[j]) - self.theoryTrue.gCl(spec,ellMids[j])) * ellWidths(j)\
                              * invcov(i,j)
        return invcov 
                
                
