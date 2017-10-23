from Utilities import *

class OneBodyMatrix:
   def __init__(self,paraDef, strOneBody, sitesDef, superLatticeVectors, mu, isNambu=False):
      #self.index1 = []
      #self.index2 = []
      #self.values = []
      #self.Rx = []
      #self.Ry = []
      #self.Rz = []
      
      sparse_0 = {} # k indep, only triangular
      sparse_k = {} # k dep: we stock Folding vector (to calculate phases), only triangular
      
      list1 = strOneBody.split('\n')
      
      for i in range(0,sitesDef.n_label): #for mu, the chemical potential is negative in the t-matrix (since we do -t, at the end it will be +mu in the green function)
         addKeyElementSparseMatrix(sparse_0,(i,i),-mu)
      
      #for jj in range(0,len(paraDef.parametersList)):
      for parameterName in paraDef.parametersDict:
         #parameterName  = paraDef.parametersList[jj].name
         #parameterValue = paraDef.parametersList[jj].value
         parameterValue = paraDef.parametersDict[parameterName]
         
         for ii in range(0,len(list1)):
            rule = list1[ii].split()
            if rule[0] == parameterName:
               ruleName, ruleDelta, ruleCoeff, ruleOption, isDelta, band1, band2 = AnalyseRule2(rule)
             
               if ruleOption == '': #basic hopping
                  assert isDelta
                  
                  for i in range(0,sitesDef.n_label):
                     if sitesDef.band[i] == band1:
                        newBand = band2
                        newPos = sitesDef.positions[i]+ruleDelta;
                        newSpin = sitesDef.spin[i]

                        foldingVector,newPos1 = Fold(array(newPos),sitesDef.sitesPositions,superLatticeVectors)   
                        newHash = sitesDef.HashSpinBandPosition(newSpin,newBand,newPos1)
                        newIndex = sitesDef.hashableSites.index(newHash)
                        if all(foldingVector == array([0,0,0])):
                           if newIndex > i:
                              addKeyElementSparseMatrix(sparse_0,(i,newIndex),parameterValue*ruleCoeff)
                           else:
                              addKeyElementSparseMatrix(sparse_0,(newIndex,i),parameterValue*ruleCoeff)
                        else:
                           if newIndex > i:
                              addKeyListSparseMatrix(sparse_k,(i,newIndex),[parameterValue*ruleCoeff,foldingVector])
                           elif newIndex < i:
                              addKeyListSparseMatrix(sparse_k,(newIndex,i),[parameterValue*ruleCoeff,-foldingVector])  
                           else:
                              addKeyListSparseMatrix(sparse_k,(i,newIndex),[parameterValue*ruleCoeff,foldingVector])
                              addKeyListSparseMatrix(sparse_k,(newIndex,i),[parameterValue*ruleCoeff,-foldingVector])  
                  
               elif ruleOption == 'diagonal': #site chemical potential
                  assert not isDelta
                  position1 = ruleDelta
                  for i in range(0,sitesDef.n_label):
                     if all(sitesDef.positions[i] == position1):
                        addKeyElementSparseMatrix(sparse_0,(i,i),parameterValue*ruleCoeff)       
               elif ruleOption == 'diagonal+': #site chemical potential
                  assert not isDelta
                  position1 = ruleDelta
                  for i in range(0,sitesDef.n_label):
                     if all(sitesDef.positions[i] == position1) and i < sitesDef.n_sites:
                        addKeyElementSparseMatrix(sparse_0,(i,i),parameterValue*ruleCoeff)
               elif ruleOption == 'diagonal-': #site chemical potential
                  assert not isDelta
                  position1 = ruleDelta
                  for i in range(0,sitesDef.n_label):
                     if all(sitesDef.positions[i] == position1) and not (i < sitesDef.n_sites):
                        addKeyElementSparseMatrix(sparse_0,(i,i),parameterValue*ruleCoeff)       
               else:
                  assert False, 'option "' + ruleOption + '" not implemented yet'
      
      if isNambu:
         for keys in sparse_0:
            i = keys[0]
            j = keys[1]
            n = sitesDef.n_sites
            assert ((i<n) and (j<n)) or ((i>=n) and (j>=n)), 'no spin flip is allowed when using option "isNambu" or "anomalous"'
            if i >= n:
               sparse_0[keys] = -sparse_0[keys]
               
         for keys in sparse_k:
            i = keys[0]
            j = keys[1]
            n = sitesDef.n_sites
            assert ((i<n) and (j<n)) or ((i>=n) and (j>=n)), 'no spin flip is allowed when using option "isNambu" or "anomalous"'
            if i >= n:
               for nn in range(0,len(sparse_k[keys])):
                  
                  sparse_k[keys][nn][0] = -sparse_k[keys][nn][0]
                     
      self.sparse_0 = sparse_0
      self.sparse_k = sparse_k  
      self.mu = mu
      self.n_label = sitesDef.n_label
      #self.nParam = len(self.Rx) 
      
      #precalculate the k-indep matrix_0
      matrix_0 =  zeros((self.n_label,self.n_label),dtype=complex)
      for keys in self.sparse_0:
         matrix_0[keys[0],keys[1]]=self.sparse_0[keys]
         matrix_0[keys[1],keys[0]]=self.sparse_0[keys]
      self.matrix_0 = matrix_0
      
      sparseK = deepcopy(self.sparse_k)
      sparse0 = deepcopy(self.sparse_0)
      
      for keys in self.sparse_k:
         if keys[0] is not keys[1]:
            newKeys = (keys[1],keys[0])
            assert newKeys not in sparseK
            sparseK[newKeys] = deepcopy(sparseK[keys])
            for nn in range(0,len(sparseK[newKeys])):
               sparseK[newKeys][nn][1] = -sparseK[newKeys][nn][1]
      
      for keys in self.sparse_0:
         if keys[0] is not keys[1]:
            newKeys = (keys[1],keys[0])
            assert newKeys not in sparse0
            sparse0[newKeys] = deepcopy(sparse0[keys])
      
      self.sparseK = sparseK    
      self.sparse0 = sparse0    
      self.hybridizationFirstMoment = self.returnHybridizationFirstMoment()

      #self.returnHybridizationSecondMoment()
      #PrettyPrintComplexMatrix(self.hybridizationFirstMoment)

#      exit()
      #print matrix_0
      #print sparse_k
      #exit()
      
      
   def returnMatrix_k(self,k_vector):
      matrix_k =  zeros((self.n_label,self.n_label),dtype=complex)
      #print self.sparse_k
      for keys in self.sparse_k:
         #print self.sparse_k[keys]
         for nn in range(0,len(self.sparse_k[keys])):
            matrix_k[keys[0],keys[1]]+=((self.sparse_k[keys])[nn][0])*exp(1j*dot(k_vector,(self.sparse_k[keys])[nn][1]))

         matrix_k[keys[1],keys[0]]=conj(matrix_k[keys[0],keys[1]]) # matrix is hermitic
      set_printoptions(precision=3)
      resultingMatrix = self.matrix_0+matrix_k
      return resultingMatrix
      ###print 'Re\n', real(resultingMatrix)
      ###print 'Im\n', imag(resultingMatrix)
      ###print self.n_label
   
   def returnHybridizationFirstMoment(self):
      
      #the goal here is to calculate the first moment: a = sum_brillouin (delta_t)^2, 
      #found from eq. (30) from Maier, Jarrel et al... 2005
      
      #we don't need to keep exponential terms, so we do the multiplication
      #and as soon as two argument of multiplying exponential doesn't cancel, 
      #the contribution to the first moment is null because of the integral over
      #the brillouin zone. 
      
      matrixHybFirstMoment  =  zeros((self.n_label,self.n_label),dtype=float)
      #matrixHybFirstMoment2 =  zeros((self.n_label,self.n_label),dtype=float)
      #matrixHybFirstMoment22 =  zeros((self.n_label,self.n_label),dtype=float)
      #matrixHybFirstMoment3 =  zeros((self.n_label,self.n_label),dtype=float)
      #print 'salut', self.sparse_k

      #self.sparse_k[keys1]

      sparseK = self.sparseK   
      sparse0 = self.sparse0
                  
      for i1 in range(0,self.n_label):
         for i2 in range(0,self.n_label):
            for dum in range(0,self.n_label):
               keys1 = (i1,dum) 
               keys2 = (dum,i2)
               if (keys1 in sparseK) and (keys2 in sparseK):
                  for nn1 in range(0,len(sparseK[keys1])):
                     for nn2 in range(0,len(sparseK[keys2])):
                        if all((sparseK[keys1][nn1][1]+sparseK[keys2][nn2][1])==array([0,0,0])) :
                           matrixHybFirstMoment[i1,i2] += ((sparseK[keys1])[nn1][0])*((sparseK[keys2])[nn2][0])
      
      return matrixHybFirstMoment
      

   def returnHybridizationSecondMoment(self):

      sparseK = self.sparseK   
      sparse0 = self.sparse0
      
      #mu = 3.1      
      matrixSelf = zeros((self.n_label,self.n_label),dtype=float)
      #matrixSelf[0,0] = -40.4
      #matrixSelf[1,2] = 1.2
      
      matrixE = zeros((self.n_label,self.n_label),dtype=float)
      
      matrixE -= matrixSelf
      
      for keys in sparse0:
         matrixE[keys[0],keys[1]] -=  sparse0[keys]

      PrettyPrintComplexMatrix(matrixE)
      #exit()
      
      matrixHybSecondMoment = zeros((self.n_label,self.n_label),dtype=float)
      sparseK = self.sparseK   
      sparse0 = self.sparse0

      for i1 in range(0,self.n_label):
         for i2 in range(0,self.n_label):
            for dum1 in range(0,self.n_label):
               for dum2 in range(0,self.n_label):
                  keys1 = (i1,dum1)
                  keys2 = (dum1,dum2)
                  keys3 = (dum2,i2)
                  if (keys1 in sparseK) and (keys2 in sparseK) and (keys3 in sparseK):
                     for nn1 in range(0,len(sparseK[keys1])):
                        for nn2 in range(0,len(sparseK[keys2])):
                           for nn3 in range(0,len(sparseK[keys3])):
                              if all((sparseK[keys1][nn1][1]+sparseK[keys2][nn2][1]+sparseK[keys3][nn3][1])==array([0,0,0])):
                                 matrixHybSecondMoment[i1,i2] += ((sparseK[keys1])[nn1][0])*((sparseK[keys2])[nn2][0])*((sparseK[keys3])[nn3][0])
                  if (keys1 in sparseK) and (keys3 in sparseK):
                     for nn1 in range(0,len(sparseK[keys1])):
                        for nn3 in range(0,len(sparseK[keys3])):
                           if all((sparseK[keys1][nn1][1]+sparseK[keys3][nn3][1])==array([0,0,0])):
                              matrixHybSecondMoment[i1,i2] -= ((sparseK[keys1])[nn1][0])*matrixE[keys2[0],keys2[1]]*((sparseK[keys3])[nn3][0])

      PrettyPrintComplexMatrix(matrixHybSecondMoment)
      
      exit()
      #print matrixHybFirstMoment
      #exit()
      return matrixHybFirstMoment
                                                

