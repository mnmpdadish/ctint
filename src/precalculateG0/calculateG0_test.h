#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "calculateG0.h"

int test_Plaquettes2x2(int verbose) {
  int Nerror=0;
  char fileName[]="testInputFiles/plaquette2x2.model";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading model from %s:\n", fileName);
  
  Model model;
  read_Model(file,&model);
  
  cMatrixFunction g0_matsubara;
  dMatrixFunction g0_tau;
  init_cMatrixFunction(&g0_matsubara, &model);
  init_dMatrixFunction(&g0_tau, &model);
  
  calculate_G0_matsubara(&g0_matsubara, &model);
  calculate_G0_tau(&g0_matsubara,&g0_tau);
  
  FILE *fileOut = fopen("green0.dat","w");
  writeToFile_dMatrixFunction(fileOut, &g0_tau, &model);
  fclose(fileOut);  
  
  FILE *fileOut2= fopen("greenW0.dat","w");
  writeToFile_cMatrixFunction(fileOut2, &g0_matsubara, &model);
  fclose(fileOut2);  
  //writeToFile_cMatrixFunction(FILE *fileOut, cMatrixFunction * cMatFun, Model * model) {
  if(!doubleEqual(g0_tau.matrices[5].data[0], -0.74174358)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[1],  0.13809594)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[3],  0.19526082)) Nerror+=1;
  //print_dMatrix(&g0_tau.matrices[5]);
  
  free_dMatrixFunction(&g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  free_Model(&model);
  
  FILE *fileIn = fopen("greenW0.dat","r");
  readFile_cMatrixFunction(fileIn, &g0_matsubara, &model);
  fclose(fileIn);  
  
  return Nerror;
}

