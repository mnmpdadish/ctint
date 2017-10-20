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



double fabs(double);
  
#define DATA_BUFFER_SIZE_1 100
#define DATA_BUFFER_SIZE_2 10000 // must be square of DATA_BUFFER_SIZE_1


int doubleEqual(double const a, double const b) {
    return fabs(a - b) < 0.000000001;
}

int complexEqual(double complex const a, double complex const b) {
    return cabs(a - b) < 0.000000001;
}

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

