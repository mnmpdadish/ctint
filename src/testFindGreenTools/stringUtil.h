//#pragma once

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


int strBeginWithToken(char *str, char *token){
  char tmpstr1[50]; memset(tmpstr1,0,50);
  sscanf(str, "%s\n", tmpstr1);
  if (strcmp(tmpstr1,token)==0) return 1;
  else return 0;
}

int countElementInStr(char const * const str, char const*const delimiters){
  //char delim1 = ' ', delim2 = '\t', delim3 = '\n';
  size_t llen = strlen(str);
  size_t dlen = strlen(delimiters);
  //printf("llen=%zu, dlen=%zu\n",llen,dlen);
  
  unsigned int i, lastDelim=-1, N=0;
  for(i=0;i<llen+1;i++){
    int j, condition = (i==llen);
    for(j=0;j<dlen;j++) condition = condition || (str[i]==delimiters[j]); //check if str[i] is any of the delilimiter characters
    if( condition ) {
      if(i-lastDelim >1) {
        if(str[lastDelim+1]=='#') return N;
        else N++;
      }
      lastDelim=i;
    }
  }  
  return N;
}


// A str should be a string of the form "1 2 3 4 5  8 3", where each number 
// should be under 256 characters
// This function will attempt to read every integer. 
// It will stop if it meet a '#' character, useful for comments at the end of line.
int readIntInStr(char const*const str, int * arrayInt, char const*const delimiters){
  //char delim1 = ' ', delim2 = '\t', delim3 = '\n';
  
  size_t llen = strlen(str);
  size_t dlen = strlen(delimiters);
  //printf("dlen=%d\n",dlen);
  
  unsigned int i, lastDelim=-1, N=0;
  char token[256];
  for(i=0;i<llen+1;i++){
    int j, condition = (i==llen);
    for(j=0;j<dlen;j++) condition = condition || (str[i]==delimiters[j]); //check if str[i] is any of the delilimiter characters
    if( condition ) {
      //printf("%d %d \n",i,llen); fflush(stdout);
      if(i-lastDelim >1) {
        if(str[lastDelim+1]=='#') return N;
        else {
        //printf("%d %d \n",i,llen); fflush(stdout);
          memset(token,0,256);
          strncpy(&token[0],&str[lastDelim+1],(i-lastDelim-1));
          //printf("%s -- %s",token, str);
          arrayInt[N]=my_atoi(token);
          N++;
        }
      }
      lastDelim=i;
    }
  }  
  return N;
}


// A parenthesisTuple should be a string of the form "(1,2,3)".
// This function test this and read the 3 integers
/*int readIntInParenthesis(char const*const parenthesisTuple, int arrayInt[3]) {
  size_t llen = strlen(parenthesisTuple);
  if(parenthesis[0]!='(' || parenthesis[0]!=')' );
      //printf("could not translate '%s' to an integer.\n",str); 
      //exit(1);
  readIntInOneLine(parenthesisTuple, &arrayInt[1], "(,)\n");
  return 0;
}*/


