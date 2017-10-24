#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "oneBodyMatrix.h"

int test_oneBodyMatrix(int verbose) {
  int Nerror=0;
  
  char fileName[]="testInputFiles/plaquette2x2.in";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading one-body operators from %s:\n", fileName);

  tMatrix tMat;
  init_tMatrix(&tMat);
  readOperators_tMatrix(file, &tMat);
  //print_tMatrix(&tMat);
  
  MultiplePositions sites;
  init_MultiplePositions(&sites);
  readSites(file, &sites, "sites");
  print_MultiplePositions(&sites,"sites");
  
  MultiplePositions superlattice;
  init_MultiplePositions(&superlattice);
  readSites(file, &superlattice, "superlattice");
  print_MultiplePositions(&superlattice,"superlattice");
  assert(superlattice.n == 3);
  
  IntPosition a;
  a.x=-11; a.y=-21; a.z=0;
  Folding folding = Fold(a, sites, superlattice);
  printf("foldingR = (%d,%d,%d)\n", folding.R.x, folding.R.y, folding.R.z);
  printf("site = (%d,%d,%d)\n",     folding.r.x, folding.r.y, folding.r.z);
  if(folding.R.x !=-12 || folding.R.y !=-22 || folding.R.z != 0) Nerror++;
  if(folding.r.x !=1   || folding.r.y !=1   || folding.r.z != 0) Nerror++;
  
  
  defineSparse_tMatrix(&tMat, &sites, &superlattice);
  print_tMatrix(&tMat);
  
  cMatrix tMatrixK, Sol1;
  init_cMatrix(&tMatrixK,sites.n);
  init_cMatrix(&Sol1,4);
  calculate_tMatrixK_2D(&tMat, &tMatrixK, 0.0, 0.25*M_PI);
  //printf("\ntMatrixK=\n"); print_cMatrix(&tMatrixK);
  
  double complex * p = Sol1.data;
  *p++=-0.8+0.0*I; *p++=-2.0+0.0*I;  *p++=-1.0-1.0*I; *p++= 0.6+0.6*I; 
  *p++=-2.0+0.0*I; *p++=-0.8+0.0*I;  *p++= 0.6+0.6*I; *p++=-1.0-1.0*I; 
  *p++=-1.0+1.0*I; *p++= 0.6-0.6*I;  *p++=-0.8+0.0*I; *p++=-2.0+0.0*I; 
  *p++= 0.6-0.6*I; *p++=-1.0+1.0*I;  *p++=-2.0+0.0*I; *p++=-0.8+0.0*I; 
  transpose_cMatrix(&Sol1); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nsample of tMatrixK(0,pi/4)=\n"); print_cMatrix(&tMatrixK);
    printf("%s\n",areEqual_cMatrix(&tMatrixK,&Sol1)? "tMatrixK==Sol1": "tMatrixK!=Sol1");
  }
  if(!areEqual_cMatrix(&tMatrixK,&Sol1)) Nerror++;

  dMatrix hybFM, Sol2;
  init_dMatrix(&hybFM,sites.n);
  init_dMatrix(&Sol2,4);
  calculate_hybFirstMoments(&tMat, &hybFM);
  
  double * q = Sol2.data;
  *q++=2.43;  *q++=-0.2;  *q++=-0.2;  *q++=-0.24; 
  *q++=-0.2;  *q++=2.43;  *q++=-0.24; *q++=-0.2; 
  *q++=-0.2;  *q++=-0.24; *q++=2.43;  *q++=-0.2; 
  *q++=-0.24; *q++=-0.2;  *q++=-0.2;  *q++=2.43; 
  transpose_dMatrix(&Sol2); // with this meth of input, we need to transpose
  
  if(verbose){
    printf("\nhybFM=\n"); print_dMatrix(&hybFM);
    //printf("\nSol2=\n");  print_dMatrix(&Sol2);
    printf("%s\n",areEqual_dMatrix(&hybFM,&Sol2)? "hybFM==Sol2": "hybFM!=Sol2");
  }
  if(!areEqual_dMatrix(&hybFM,&Sol2)) Nerror++;
  
  free_cMatrix(&tMatrixK);
  free_dMatrix(&hybFM);
  free_tMatrix(&tMat);
  free_MultiplePositions(&sites);
  free_MultiplePositions(&superlattice);

  
  return Nerror;
}


