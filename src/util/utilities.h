#pragma once

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

#include "stringUtil.h"


/*// Constants rounded for 21 decimals. 
#define M_E 2.71828182845904523536
#define M_LOG2E 1.44269504088896340736
#define M_LOG10E 0.434294481903251827651
#define M_LN2 0.693147180559945309417
#define M_LN10 2.30258509299404568402
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define M_PI_4 0.785398163397448309616
#define M_1_PI 0.318309886183790671538
#define M_2_PI 0.636619772367581343076
#define M_1_SQRTPI 0.564189583547756286948
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2 1.41421356237309504880
#define M_SQRT_2 0.707106781186547524401
*/

double fabs(double);
  
//#define DATA_BUFFER_SIZE_1 100
//#define DATA_BUFFER_SIZE_2 10000 // must be square of DATA_BUFFER_SIZE_1
#define INIT_CAPACITY 16

int doubleEqual(double const a, double const b) {
    return fabs(a - b) < 0.000001;
}

int complexEqual(double complex const a, double complex const b) {
    return cabs(a - b) < 0.000001;
}

void readDouble(FILE * file, char * name,  double * value) {

    rewind(file);
    
    char tempbuff[200];
    while(!feof(file)) 
    {
        if (fgets(tempbuff,200,file))
        {   
            int lenString = strlen(tempbuff);
            if(tempbuff[lenString-1]!='\n') {
              printf("error. buffer too small. %c\n",tempbuff[lenString]);
              exit(1);
            }
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atof(tmpstr2); return;} 
        }
    }
    printf("\ncannot find the %s parameter.\n", name);
    exit(1);
}

void readInt(FILE * file, char * name,  int * value) {

    rewind(file);
    
    char tempbuff[200];
    while(!feof(file)) 
    {
        if (fgets(tempbuff,200,file))
        {
            int lenString = strlen(tempbuff);
            if(tempbuff[lenString-1]!='\n') {
              printf("error. buffer too small. %c\n",tempbuff[lenString]);
              exit(1);
            }
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atoi(tmpstr2); return;} 
        }
    }
    printf("\ncannot find the %s parameter.\n", name);
    exit(1);
}




int countLineFlag(FILE * file, char *flag) {

  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, N=0;
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        if(strBeginWithToken(tempbuff,flag)) found=1; 
      }
      else{
        if(tempbuff[0] == '#') continue;
        else if((tempbuff[0] != '\n') && countElementInStr(tempbuff," \t\n") != 0) N++;
        else break;
      }
    }
    //printf("\n");
  }
  return N;
}


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

