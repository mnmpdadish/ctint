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

FILE * fopenSafe(char fileName[], char mode[], int verbose){
  FILE * file = fopen(fileName, mode);
  if(file == NULL) {
    printf("error: file %s not found.\nterminated.\n", fileName); 
    exit(1);
  }
  else {
    if(verbose) printf("opening file %s\n", fileName);
  }
  return file;
}


void copyFile(char fileNameIn[], char fileNameOut[], const char * paramName[], const char * replacementValue[], unsigned int nElement) {

  FILE * fileIn = fopenSafe(fileNameIn, "r", 1);
  FILE * fileOut= fopenSafe(fileNameOut, "w", 1);
  
  unsigned int n1,n2,n,i;
  char tempbuff[2048];
  while(!feof(fileIn)) 
  {
      if (fgets(tempbuff,2048,fileIn))
      {
          char tmpstr1[256];
          char tmpstr2[256];
          sscanf(tempbuff, "%255s%n %n%255s\n", tmpstr1, &n1, &n2, tmpstr2);
          int found=0;
          for(i=0;i<nElement;i++){
            if (strcmp(tmpstr1,paramName[i])==0) {
              found=1;
              fprintf(fileOut,"%s", paramName[i]);
              n=n2-n1;
              while(n>0){
                fprintf(fileOut," ");
                n--;
              }
              fprintf(fileOut,"%s\n", replacementValue[i]);
            }
          }
          if(!found) fprintf(fileOut,"%s", tempbuff);
      }
  }
  fclose(fileIn);
  fclose(fileOut);
}

void main() {

  char test[256];
  sprintf(test, "%f", 1.985); // puts string into buffer
  
  
  const char* paramName[] = { "mu", "measure_i", "Jim" };
  const char* replacementValue[] = { test, "67", "Jim" };
  copyFile("params0", "params1", paramName, replacementValue, 3); 
}
