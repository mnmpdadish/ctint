#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "findGreenSymmetries.h"

int test_readSymmetriesPlaquettes2x2(int verbose, char fileName[]) {
  int nSites=4;
  //char fileName[]="testInputFiles/plaquette2x2.in";
  FILE * file = fopenSafe(fileName, "rt", verbose);
  
  int nSym=countLineFlag(file, "symmetry_generators");
  if(verbose) printf("nsym=%d\n", nSym);
  
  Symmetries sym; initSymmetries(&sym, nSym, nSites); //kind of a constructor.
  readSymmetries(file, nSites, &sym, "symmetry_generators");
  if(verbose) printSymmetries(&sym);
  
  GreenSymmetriesMatrix greenSymMat;
  initGreenSymmetriesMatrix(&greenSymMat,nSites,&sym);
  
  if(verbose){
    printf("\nsymmetrized matrix:\n");
    printGreenSymmetriesMatrix(&greenSymMat);
    printf("\n");
  }
  
  int solution_i[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int solution_j[] = {0,1,1,3,1,0,3,1,1,3,0,1,3,1,1,0};
  int Nerror=0,i,j,N=greenSymMat.nSites;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if((greenSymMat.i[N*i+j]!=solution_i[N*i+j]) || (greenSymMat.j[N*i+j]!=solution_j[N*i+j])) Nerror++;
    }
  } 
  
  
  
  freeGreenSymmetriesMatrix(&greenSymMat);
  freeSymmetries(&sym);

  return Nerror;
}

int test_readSymmetriesPlaquettes4x4(int verbose, char fileName[]) {
  int nSites=16;
  //char fileName[]="testInputFiles/plaquette4x4.in";
  FILE * file = fopenSafe(fileName, "rt", verbose);
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  if(verbose) printf("\nreading symmetries from %s:\n", fileName);

  //Symmetries symmetries
  int nSym=countLineFlag(file, "symmetry_generators");
  if(verbose) printf("nsym=%d\n", nSym);
  
  Symmetries sym; initSymmetries(&sym, nSym, nSites); //kind of a constructor.
  readSymmetries(file, nSites, &sym, "symmetry_generators");
  if(verbose) printSymmetries(&sym);
  
  GreenSymmetriesMatrix greenSymMat;
  initGreenSymmetriesMatrix(&greenSymMat,nSites,&sym);
  
  if(verbose){
    printf("\nsymmetrized matrix:\n");
    printGreenSymmetriesMatrix(&greenSymMat);
  }
  
  int solution_i[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,0,
                      1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,5,5,
                      1,1,5,5,1,0,1,1,0,0,1,1,0,1,5,5,1,1,5,5,1,0,1,1,0,0,1,1,0,
                      1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,
                      1,1,0,1,5,5,1,1,5,5,1,0,1,1,0,0,1,1,0,1,5,5,1,1,5,5,1,0,1,
                      1,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                      0,0,0,0,0,0,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,
                      1,1,1,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  int solution_j[] = {0,1,2,3,1,5,6,7,2,6,10,11,3,7,11,15,1,1,2,2,4,5,6,7,7,9,
                      10,11,7,13,14,11,2,2,1,1,7,6,5,4,11,10,9,7,11,14,13,7,
                      3,2,1,0,7,6,5,1,11,10,6,2,15,11,7,3,1,4,7,7,1,5,9,13,2,6,
                      10,14,2,7,11,11,5,5,6,6,5,5,6,9,6,6,10,10,6,9,10,10,6,6,
                      5,5,9,6,5,5,10,10,6,6,10,10,9,6,7,7,4,1,13,9,5,1,14,10,
                      6,2,11,11,7,2,2,7,11,11,2,6,10,14,1,5,9,13,1,4,7,7,6,
                      9,10,10,6,6,10,10,5,5,6,9,5,5,6,6,10,10,9,6,10,10,6,6,
                      9,6,5,5,6,6,5,5,11,11,7,2,14,10,6,2,13,9,5,1,7,7,4,
                      1,3,7,11,15,2,6,10,11,1,5,6,7,0,1,2,3,7,13,14,11,7,9,10,
                      11,4,5,6,7,1,1,2,2,11,14,13,7,11,10,9,7,7,6,5,4,2,2,1,1,
                      15,11,7,3,11,10,6,2,7,6,5,1,3,2,1,0};

  int Nerror=0,i,j,N=greenSymMat.nSites;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if((greenSymMat.i[N*i+j]!=solution_i[N*i+j]) || (greenSymMat.j[N*i+j]!=solution_j[N*i+j])) Nerror++;
    }
  } 
  
  if(verbose) {
    printf("nIndep=%d\n",greenSymMat.nIndep);
    printf("\n");
  }
  
  freeGreenSymmetriesMatrix(&greenSymMat);
  freeSymmetries(&sym);

  return Nerror;
}


int test_arbitrary(int verbose) {
  int nSites=16;
  char fileName[]="testInputFiles/plaquette4x4_test.in";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading symmetries from %s:\n", fileName);

  //Symmetries symmetries
  int nSym=countLineFlag(file, "symmetries");
  printf("nsym=%d\n", nSym);
  
  Symmetries sym; initSymmetries(&sym, nSym, nSites); //kind of a constructor.
  readSymmetries(file, nSites, &sym, "symmetries");
  printSymmetries(&sym);
  
  GreenSymmetriesMatrix greenSymMat;
  initGreenSymmetriesMatrix(&greenSymMat,nSites,&sym);
  
  if(verbose){
    printf("\nsymmetrized matrix:\n");
    printGreenSymmetriesMatrix(&greenSymMat);
    printf("nIndep=%d\n",greenSymMat.nIndep);
    printf("\n");
  }
  
  freeGreenSymmetriesMatrix(&greenSymMat);
  freeSymmetries(&sym);

  return 0;
}

