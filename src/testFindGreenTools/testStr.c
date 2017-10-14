#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int countNumberOfElementInOneLine(const char *line){
  char delim1 = ' ', delim2 = '\t';
  size_t llen = strlen(line);
  
  unsigned int i, lastDelim=-1, N=0;
  for(i=0;i<llen;i++){
    if((line[i]==delim1) || (line[i]==delim2)){
      if(i>0 && i-lastDelim >1) {
        N++;
      }
      lastDelim=i;
    }
  }  
  return N;
}

int readIntInOneLine(const char *line, int * arrayInt){
  char delim1 = ' ', delim2 = '\t';
  size_t llen = strlen(line);
  
  unsigned int i, lastDelim=-1, N=0;
  char token[256];
  for(i=0;i<llen;i++){
    if((line[i]==delim1) || (line[i]==delim2)){
      if(i>0 && i-lastDelim >1) {
        memset(token,0,strlen(token));
        strncpy(&token[0],&line[lastDelim+1],(i-lastDelim-1));
        //printf("%d %d ",i, lastDelim );
        //printf("%s \n",token );
        arrayInt[N]=atoi(token);
        N++;
      }
      lastDelim=i;
    }
  }  
  return N;
}


int main (int argc, char **argv)
{
  char line[] = "1 3 2 \t  44 5 ";
  unsigned int nElement = countNumberOfElementInOneLine(line);
  int *arrayInt = (int *) malloc(nElement * sizeof(int)); 
  readIntInOneLine(line, arrayInt);
  int i;
  for(i=0;i<nElement;i++) printf("%d ", arrayInt[i]);
  return 0;
}
