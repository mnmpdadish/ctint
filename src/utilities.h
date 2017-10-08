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


// reminder:
// int * ptr;
// int const * ptrToConst;
// int * const constPtr;
// int const * const constPtrToConst;
// At first, I added const AFTER each *, but in the end, it is just more confusing than nothing.


int dgemm_(char const*, char const*, int const*, int const*, int const*,
                  double const*, double const*, int const*, double const*, 
                  int const*, double const*, double *, int const*);
int dswap_(int const*, double*, int const*, double*, int const*);
int dgetrf_(int const*, int const*, double const*, int const*,int*,int*);
int dgetri_(int const*, double*, int const*, int const*, double*, int const*, int*);
int dgemv_(char const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
double ddot_(int const*, double const*, int const*, double const*, int const*);
int dger_(int const*, int const*, double const*, double const*, int const*, double const*, int const*, double *, int const*);
int dcopy_(int const*, double const*, int const* , double*, int const*);
int daxpy_(int const*, double const*, double const*, int const*, double*, int const*);
int dscal_(int const*, double const*, double*, int const*);

double fabs(double);
  
#define dataBufferSize1 100
#define dataBufferSize2 10000 // must be square of dataBufferSize1


void readDouble(FILE * file, char * name,  double * value) {

    rewind(file);
    
    char tempbuff[200];
    while(!feof(file)) 
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atof(tmpstr2); return;} 
        }
    }
    printf("\ncannot find the %s parameter in 'model.dat'", name);
    exit(1);
}

void readInt(FILE * file, char * name,  int * value) {

    rewind(file);
    
    char tempbuff[200];
    while(!feof(file)) 
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atoi(tmpstr2); return;} 
        }
    }
    printf("\ncannot find the %s parameter in 'model.dat'", name);
    exit(1);
}

