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
  print_tMatrix(&tMat);
  free_tMatrix(&tMat);

  MultiplePositions sites;
  init_MultiplePositions(&sites);
  readSites(file, &sites, "sites");
  print_MultiplePositions(&sites,"sites");
  free_MultiplePositions(&sites);

  MultiplePositions lattice;
  init_MultiplePositions(&lattice);
  readSites(file, &lattice, "lattice");
  print_MultiplePositions(&lattice,"lattice");
  free_MultiplePositions(&lattice);

  MultiplePositions superlattice;
  init_MultiplePositions(&superlattice);
  readSites(file, &superlattice, "superlattice");
  print_MultiplePositions(&superlattice,"superlattice");
  free_MultiplePositions(&superlattice);

  
  return Nerror;
}


