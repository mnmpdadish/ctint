#pragma once

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <assert.h>
#include <ctype.h>


int my_atoi(char const * const str) {
  //char strArray[256];
  //char *str2 = &strArray[0];
  //assert(strlen(str)<256);
  //strcpy(str2, str);
  
  int i=0, res=0, minus = str[i] == '-';
  if (minus) i++;
  
  int N=strlen(str); //important to not evaluate each iteration
  while(i<N){  
    //printf("strlen:%d %d %s\n",i,N,str);
    if(!isdigit(str[i])) {
      printf("could not translate '%s' to an integer.\n",str); 
      exit(1);
    }
    res = res*10 + (str[i++] - '0');
  }

  
  return minus ? -res : res;
}

int countNumberOfElementInOneLine(char const * const line){
  char delim1 = ' ', delim2 = '\t', delim3 = '\n';
  size_t llen = strlen(line);
  
  unsigned int i, lastDelim=-1, N=0;
  for(i=0;i<llen+1;i++){
    if((i==llen) || (line[i]==delim1) || (line[i]==delim2) || (line[i]==delim3)){
      if(i>0 && i-lastDelim >1) {
        N++;
      }
      lastDelim=i;
    }
  }  
  return N;
}

int readIntInOneLine(char const*const line, int * arrayInt){
  char delim1 = ' ', delim2 = '\t', delim3 = '\n';
  size_t llen = strlen(line);
  
  unsigned int i, lastDelim=-1, N=0;
  char token[256];
  for(i=0;i<llen+1;i++){
    if( (i==llen) || (line[i]==delim1) || (line[i]==delim2) || (line[i]==delim3)){
      //printf("%d %d \n",i,llen); fflush(stdout);
      if(i-lastDelim >1) {
        //printf("%d %d \n",i,llen); fflush(stdout);
        memset(token,0,256);
        strncpy(&token[0],&line[lastDelim+1],(i-lastDelim-1));
        //printf("%s -- %s",token, line);
        arrayInt[N]=my_atoi(token);
        N++;
      }
      lastDelim=i;
    }
  }  
  return N;
}


int test_readIntInOneLine()
{
  char line[] = "1 3 2 \t \t\t  34784 5 ";
  unsigned int nElement = countNumberOfElementInOneLine(line);
  int *arrayInt = (int *) malloc(nElement * sizeof(int)); 
  readIntInOneLine(line, arrayInt);
  int i;
  for(i=0;i<nElement;i++) printf("%d ", arrayInt[i]);
  return 0;
}
