#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../util/stringUtil.h"
#include "../util/utilities.h"
#include "../util/arrays/array.h"

//#define DATA_BUFFER_SIZE_W 200

typedef struct {
  Array_int *permutations;
  unsigned int n;
  //unsigned int timeReversal;
} Symmetries;

int initSymmetries(Symmetries * sym, int nSym, int nSites) {
  //assert(nSites<=64);
  sym->permutations = (Array_int *) malloc(nSym * sizeof (Array_int ));
  sym->n=nSym;
  int i;
  for(i=0;i<nSym;i++) Array_int_init(&sym->permutations[i]);
  return 0;
}

int freeSymmetries(Symmetries * sym) {
  int i;
  for(i=0;i<sym->n;i++) Array_int_free(&sym->permutations[i]);
  free(sym->permutations);
  return 0;
}


int addElementToPermutation(Array_int *perm, unsigned int value, unsigned int nSites) {
  unsigned int i;
  //the four next lines are there to ensure that the permutation has only one of 
  //each number, between 0 and nSites-1.
  for(i=0;i<perm->size;i++){
    assert(perm->data[i] != value );
  }
  //printf("%d %d\n", value, perm->nSites);
  assert(value<nSites && value>=0);

  Array_int_push(perm, value);
  
   //perm->data[perm->size]=value;
  //perm->nAssigned++;
  //assert(perm->nAssigned <= perm->nSites);
  return 0;
}


char charHexa(unsigned int digit) {
  //any cluster should not have more than 64 sites, even in 2050 xD
  char characters[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!?"; 
  char char1;
  if(digit < 64 && digit >= 0) char1 = characters[digit];
  else char1 = '.';
  return char1;
}


int printSymmetries(Symmetries * sym){
  int i,j;
  for(i=0;i<sym->n;i++){
    for(j=0;j<sym->permutations[i].size;j++){
      printf("%c", charHexa(sym->permutations[i].data[j]));
    }
    printf("\n");
  } 
  return 0;
}
      


int readSymmetries(FILE * file, int nSites, Symmetries *sym, char *flag) {
  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, i=0;
  //char name[]="symmetries";
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      //printf("%d seen= %s",found, tempbuff);
      if(found==0){
        if(strBeginWithToken(tempbuff,flag)) found=1; 
      }
      else{
        unsigned int nElement = countElementInStr(tempbuff, " \t\n");
        if(tempbuff[0] == '#') continue;
        //else if(strBeginWithToken(tempbuff,"time")) sym->timeReversal=1;
        else if((tempbuff[0] != '\n') && nElement != 0) {
		      int arrayInt[nElement];
          readIntInStr(tempbuff, arrayInt, " \t\n");
          int j;
          for(j=0;j<nElement;j++) {
            //printf("number= %d\n",arrayInt[j]);
            addElementToPermutation(&(sym->permutations[i]),arrayInt[j],nSites);
          }
          assert(nSites == nElement);
          i++;
        }
        else break;
      }
    }
  }
  return 0;
}





typedef struct {
  unsigned int *i;
  unsigned int *j;
  unsigned int nElement;
  unsigned int nSites;
  unsigned int nIndep;
  double * greenfunctionsTAU;
  //double * greenfunctionsIWn_Re;
  //double * greenfunctionsIWn_Im;
  unsigned int capacity;
  //Array_complex *greenfunctionsTAU;
} GreenMatrix;


int printGreenMatrix(GreenMatrix * greenMatrix) {  
  int i,j,N=greenMatrix->nSites;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      printf("%c%c ",charHexa(greenMatrix->i[N*i+j]),charHexa(greenMatrix->j[N*i+j]));
    }
    printf("\n");
  } 
  return 0;
}




int countIndependantGreen(GreenMatrix * greenMatrix) {  
  int i,j,N=greenMatrix->nSites, nIndep=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(greenMatrix->i[N*i+j]==i && greenMatrix->j[N*i+j]==j) nIndep++;
    }
  } 
  return nIndep;
}



int symmetrizeOneGreenElement(GreenMatrix * greenMatrix, Symmetries * sym) {  
  int i,j,k, N=greenMatrix->nSites;
  for(k=0; k<sym->n; k++){
    for(i=0;i<greenMatrix->nSites;i++){
      for(j=0;j<greenMatrix->nSites;j++){
        int permutation_i = i;
        int permutation_j = j;
        int isFirst=1, n=0;
        while( ((permutation_i != i) || (permutation_j != j)) || isFirst) { //while didn't loop on orbit
          n++;
          permutation_i= sym->permutations[k].data[permutation_i];
          permutation_j= sym->permutations[k].data[permutation_j];
          int refIndex=N*i+j ;
          int newIndex=N*permutation_i+permutation_j;
          if( (greenMatrix->i[newIndex] != greenMatrix->i[refIndex]) || (greenMatrix->j[newIndex] != greenMatrix->j[refIndex]) ) {
            //printf("k=%d  i=%d j=%d %d %d  %c%c %c%c  %d %d\n",k,i,j,permutation_i,permutation_j, charHexa(greenMatrix->i[refIndex]), charHexa(greenMatrix->j[refIndex]), charHexa(greenMatrix->i[newIndex]), charHexa(greenMatrix->j[newIndex]), refIndex, newIndex);
            if( N*greenMatrix->i[refIndex] + greenMatrix->j[refIndex] < N*greenMatrix->i[newIndex] + greenMatrix->j[newIndex] ) {
            //refIndex < newIndex){
              greenMatrix->i[newIndex]=greenMatrix->i[refIndex];
              greenMatrix->j[newIndex]=greenMatrix->j[refIndex];
            }
            else{
              greenMatrix->i[refIndex]=greenMatrix->i[newIndex];
              greenMatrix->j[refIndex]=greenMatrix->j[newIndex];
            }
            //printf("k=%d  i=%d j=%d %d %d  %c%c %c%c\n",k,i,j,permutation_i,permutation_j, charHexa(greenMatrix->i[refIndex]), charHexa(greenMatrix->j[refIndex]), charHexa(greenMatrix->i[newIndex]), charHexa(greenMatrix->j[newIndex]));

            return 1; //every change result in a function exit with 1
          }
          isFirst=0;
          if(n>20) exit(1);
        }
      }
    }
  }
  return 0; //once there is no more changes, this function should exit with 0
}


int initGreenMatrix(GreenMatrix * greenMatrix, int nSites, Symmetries *sym) {

  
  greenMatrix->i = (unsigned int *) malloc(nSites*nSites * sizeof (unsigned int));
  greenMatrix->j = (unsigned int *) malloc(nSites*nSites * sizeof (unsigned int));
  greenMatrix->nElement = nSites*nSites;
  greenMatrix->nSites = nSites;
  
  int i,j;
  for(i=0;i<nSites;i++){
    for(j=0;j<nSites;j++){
      greenMatrix->i[i*nSites+j]=(i<j)? i:j;
      greenMatrix->j[i*nSites+j]=(i<j)? j:i;
    }
  } 
  
  //if(verbose){
  //  printf("\nmost General matrix:\n");
  //  printGreenMatrix(greenMatrix);
  //}
  
  while(symmetrizeOneGreenElement(greenMatrix,sym));  //this line symmetrize the matrix
  
  
  greenMatrix->nIndep = countIndependantGreen(greenMatrix);
  
  //if(verbose) printf("nIndep=%d\n",nGreen);
  return 0;
}

int freeGreenMatrix(GreenMatrix * greenMatrix) {
  free(greenMatrix->i);
  free(greenMatrix->j);
  return 0;
}


/*
int main() {
  int nSites=16;
  char fileName[]="plaquette4x4.in";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found", fileName); exit(1);}
  printf("reading symmetries from %s:\n", fileName);

  //Symmetries symmetries
  int nSym=countLineFlag(file, "symmetries");
  printf("nsym=%d\n", nSym);
  
  Symmetries sym; buildSymmetries(&sym, nSym, nSites); //kind of a constructor.
  readSymmetries(file, nSites, &sym);
  printSymmetries(&sym);
  
  GreenMatrix greenMatrix;
  buildGreenMatrix(&greenMatrix,nSites,&sym);
  
  // do stuff with green.
  
  freeGreenMatrix(&greenMatrix);
  freeSymmetries(&sym);

  return 0;
}
*/

