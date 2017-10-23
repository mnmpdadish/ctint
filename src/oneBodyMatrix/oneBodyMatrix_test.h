#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "oneBodyMatrix.h"

int test_oneBodyMatrix(int verbose) {
  int Nerror=1;
  
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
  
  /*
  MultiplePositions lattice;
  init_MultiplePositions(&lattice);
  readSites(file, &lattice, "lattice");
  print_MultiplePositions(&lattice,"lattice");
  assert(lattice.n > 0 && lattice.n <= 3);
  */
  
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
    
  free_tMatrix(&tMat);
  free_MultiplePositions(&sites);
  //free_MultiplePositions(&lattice);
  free_MultiplePositions(&superlattice);

  
  return Nerror;
}


