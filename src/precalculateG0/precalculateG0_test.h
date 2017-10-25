#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "precalculateG0.h"

int test_Plaquettes2x2(int verbose) {
  int Nerror=1;
  char fileName[]="testInputFiles/plaquette2x2.model";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading model from %s:\n", fileName);
  
  Model model;
  read_Model(file,&model);
  
  IndepFunctionComplex g0;
  init_IndepFunctionComplex(&g0, &model);
  precalculate_G0(&g0, &model);
  
  FILE *fileOut = fopen("green0.dat","w");
  writeToFile_IndepFunctionComplex(fileOut, &g0);
  
  free_IndepFunctionComplex(&g0);
  free_Model(&model);
  
  return Nerror;
}

