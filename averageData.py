#!/usr/bin/python


import os,sys
#import numpy as np

def print_usage():
  print('\n wrong number of input parameter.\n\nexample of usage:\n$ ./average.py self 5\n\nto average the self*.dat file for iteration 5 and above (stop when does not find the file).\n')
  exit()

def ReadFile(fileName):
  f = open(fileName, 'r')
  file1 = f.readlines()
  f.close
  
  file2 = []
  for lines in file1:
    file2.append(lines.strip())
  return file2
   
input_string_len = len(sys.argv[:])

if(input_string_len !=3): 
  print(sys.argv[2])
  print_usage()


print(sys.argv[:])
list_of_fileList = []
nFile=input_string_len-1
nLine = 0

fileExists = True
iteration = int(sys.argv[2])
while(fileExists):
  
##for ii in range(1,input_string_len):
  fileName = './'+sys.argv[1]+str(iteration)+'.dat'
  if os.path.isfile(fileName):
    print fileName
    fileList = ReadFile(fileName)
    list_of_fileList.append(fileList)
    if nLine ==0:
      nLine = len(fileList)
    else:
      if nLine != len(fileList):
        print("error: files must be of same nature and same number of line.\nterminated.\n")
        exit()
    iteration+=1
  else:
    fileExists=False


#nElement = 0
nElement = len(list_of_fileList[0][0].split()) - 1 #minus one because of the # character
if(nElement % 2 != 1): #must be odd
  print('wrong number of columns')
  exit()



fileName = './'+sys.argv[1]+'_averaged_'+ sys.argv[2]+'-'+str(iteration-1)+'.dat'
print('\noutput: '+fileName+'\n')
fileOut = open(fileName,'w')
for ii in range(nLine):
  lineArray = [0.0]*nElement
  if list_of_fileList[0][ii][0] == '#':
    #print list_of_fileList[0][ii]
    line = list_of_fileList[0][ii] + '\n'
  else:
    for jj in range(nFile):
      #print(list_of_fileList[jj][ii][0])
      lineList = list_of_fileList[jj][ii].split()
      for kk in range(nElement):
        lineArray[kk] += float(lineList[kk])/nFile
        
    line = "% 3.6e  " % lineArray[0]
    for kk in range(nElement-2):
      line +="% 3.6e % 3.6e  " % (lineArray[kk+1],lineArray[kk+2])
    line += "\n"
  #print(lineArray)
  #print(line)
  fileOut.write(line)

fileOut.close()

exit()

