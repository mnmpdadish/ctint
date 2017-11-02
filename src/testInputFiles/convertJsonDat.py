import json
import os
from math import *
from numpy import *

def ReadJson(fileName):
   fp = open(fileName, 'r')
   jsonDict = json.load(fp)
   fp.close()
   return jsonDict
   
def writeDataFile(realVectorX, complexVectorYList,fileName):
   if os.path.isfile(fileName):
      print 'warning, overwriting %s' % fileName
   fp = open(fileName, 'w')
   #fp.write(header)  #header
   isCompl = True #any(iscomplex(complexVectorYList))
   for ii in range(0,len(realVectorX)):
      lineToPrint = '%  10.7f' % realVectorX[ii]
      lineToPrint += ' % 9.7f' % real(complexVectorYList[ii])
      lineToPrint += ' % 9.7f' % imag(complexVectorYList[ii])         
      fp.write(lineToPrint + '\n')
   fp.close()
   

def main1():
   hyb = ReadJson('Hyb8.json')
   #################################################
   #WRITE YOUR CODE HERE

   omega_n=[]
   complexVectorYList=[]
   for nn in range(0,len(hyb['A']['imag'])):
      omega_n.append((2*nn+1)*3.1416/hyb['A']['beta'])
      complexVectorYList.append(hyb['A']['real'][nn] + 1.0j* hyb['A']['imag'][nn])
      
   print omega_n
   print complexVectorYList
   #################################################
   writeDataFile(omega_n, complexVectorYList,'hyb8.dat')
   
main1()
