#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../util/stringUtil.h"
#include "../util/arrays/array.h"

typedef struct {
  char label[64]; //name should be under 64 characters.
  double coefficient;
  int posDiff3D[3]; 
} OperatorDef;


typedef struct {
  unsigned int n;
  unsigned int capacity;
  OperatorDef * operators;
} tMatrix;


int init_tMatrix(tMatrix * tMat) {
  tMat->capacity = INIT_CAPACITY;
  tMat->operators = (OperatorDef *) malloc(tMat->capacity * sizeof (OperatorDef ));
  tMat->n=0;
  return 0;
}

void addOperators_tMatrix(tMatrix * tMat, char label[64], double coefficient, int posDiff3D[3]) {
  if(tMat->n > tMat->capacity) tMat->operators = realloc(tMat->operators, (tMat->capacity *= 2) * sizeof(OperatorDef));
  memset(tMat->operators[tMat->n].label,0,64);
  strcpy(tMat->operators[tMat->n].label,label);
  tMat->operators[tMat->n].coefficient=coefficient;
  unsigned int i;
  for(i=0;i<3;i++) tMat->operators[tMat->n].posDiff3D[i]=posDiff3D[i];
  tMat->n++;  
}

  

void readOperators_tMatrix(FILE * file, tMatrix * tMat){
  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0;
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        if(strBeginWithToken(tempbuff,"one-body")) found=1; 
      }
      else{
        unsigned int nElement = countElementInStr(tempbuff, " \t\n");
        if(tempbuff[0] == '#') continue;
        else if(((tempbuff[0] != '*') || (tempbuff[0] != '\n')) && nElement != 0) {
          
          char label[64];
          float coefficient;
          int posDiff3D[3];
          int nRead;
          nRead = sscanf(tempbuff,"%s (%d,%d,%d) %f \n",label, &posDiff3D[0],&posDiff3D[1],&posDiff3D[2],&coefficient);
          addOperators_tMatrix(tMat, label, coefficient, posDiff3D);
          
          if(nRead!=5) {
            printf("Cannot read correctly the one-body line: \n%s", tempbuff); 
            exit(1);
          }
        }
        else break;
      }
    }
  }
}



void print_tMatrix(tMatrix * tMat) {
  unsigned int i;
  for(i=0;i<tMat->n;i++) {
    printf("operator %d: '%s' with coeff %f and position difference vector: (%d,%d,%d)\n", 
            i, tMat->operators[i].label, tMat->operators[i].coefficient,
            tMat->operators[i].posDiff3D[0],
            tMat->operators[i].posDiff3D[1],
            tMat->operators[i].posDiff3D[2]);
  }
}


void free_tMatrix(tMatrix * tMat) {
  free(tMat->operators);
}


typedef struct {
  int x;
  int y;
  int z;
} IntPosition;

typedef struct {
  unsigned int n;
  unsigned int capacity;
  IntPosition * positions;
} MultiplePositions;


int init_MultiplePositions(MultiplePositions * mPositions) {
  mPositions->capacity = INIT_CAPACITY;
  mPositions->positions = malloc(mPositions->capacity * sizeof (IntPosition));
  mPositions->n=0;
  return 0;
}

void addSites_MultiplePositions(MultiplePositions * mPositions, int pos3D[3]) {
  unsigned int n=mPositions->n;
  if(mPositions->n > mPositions->capacity) 
    mPositions->positions = realloc(mPositions->positions, (mPositions->capacity *= 2) * sizeof(IntPosition));
  mPositions->positions[n].x=pos3D[0];
  mPositions->positions[n].y=pos3D[1];
  mPositions->positions[n].z=pos3D[2];
  mPositions->n++;  
}

  

void readSites(FILE * file, MultiplePositions * mPositions, char * flag){
  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0;
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        if(strBeginWithToken(tempbuff,flag)) found=1; 
      }
      else{
        unsigned int nElement = countElementInStr(tempbuff, " \t\n");
        if(tempbuff[0] == '#') continue;
        else if(((tempbuff[0] != '*') || (tempbuff[0] != '\n')) && nElement != 0) {
          
          int pos3D[3];
          int nRead;
          nRead = sscanf(tempbuff,"(%d,%d,%d) \n", &pos3D[0],&pos3D[1],&pos3D[2]);
          addSites_MultiplePositions(mPositions, pos3D);
          
          if(nRead!=3) {
            printf("Cannot read correctly the %s line: \n%s", flag, tempbuff); 
            exit(1);
          }
        }
        else break;
      }
    }
  }
}



void print_MultiplePositions(MultiplePositions * mPositions, char *flag) {
  unsigned int i;
  for(i=0;i<mPositions->n;i++) {
    printf("%s %d with position: (%d,%d,%d)\n", flag, i, mPositions->positions[i].x
                                                       , mPositions->positions[i].y
                                                       , mPositions->positions[i].z);
  }
}


void free_MultiplePositions(MultiplePositions * mPositions) {
  free(mPositions->positions);
}

