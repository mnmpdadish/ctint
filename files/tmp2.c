#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>

void fappend(char fileName[], char header[], int iteration, double value, int verbose){
  short exist=0; 
  if( access( fileName, F_OK ) != -1 ) exist =1; 
  
  FILE * file = fopen(fileName, "a");
  if(file == NULL) {
    printf("error: cannot opent file %s.\nterminated.\n", fileName); 
    exit(1);
  }
  else {
    if(verbose) printf("opening file %s\n", fileName);
  }
  if(!exist) fprintf(file, "%s\n",header);
  fprintf(file, "%2d % 7.6f\n",iteration,value);
  fclose(file);
}



void main() {
  unsigned int i;
  for(i=0;i<200;i++) fappend("salut.txt", "#comment ca va", i, sqrt(((double) i)), 0);
  printf("salut\n");
}
