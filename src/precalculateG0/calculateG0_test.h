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
  cMatrixFunction g0_matsubara_Read;
  dMatrixFunction g0_tau;
  init_cMatrixFunction(&g0_matsubara, &model);
  init_cMatrixFunction(&g0_matsubara_Read, &model);
  init_dMatrixFunction(&g0_tau, &model);
  
  calculate_G0_matsubara(&g0_matsubara, &model, NULL);
  FILE *fileOut2= fopen("greenW0.dat","w");
  if(fileOut2 == NULL) {printf("file %s not found\n", "green0.dat"); exit(1);}
  writeToFile_cMatrixFunction(fileOut2, &g0_matsubara, &model);
  fclose(fileOut2);  


  calculate_G0_tau(&g0_matsubara,&g0_tau);
  
  FILE *fileOut = fopen("green0.dat","w");
  if(fileOut == NULL) {printf("file %s not found\n", "green0.dat"); exit(1);}
  writeToFile_dMatrixFunction(fileOut, &g0_tau, &model);
  fclose(fileOut);  
  
  
  if(!doubleEqual(g0_tau.matrices[5].data[0], -0.74174358)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[1],  0.13809594)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[3],  0.19526082)) Nerror+=1;
  //print_dMatrix(&g0_tau.matrices[5]);
  
  
  FILE *fileIn = fopen("greenW0.dat","r");
  if(fileIn == NULL) {printf("file %s not found\n", "greenW0.dat"); exit(1);}
  readFile_cMatrixFunction(fileIn, &g0_matsubara_Read, &model);
  fclose(fileIn);  
  
  FILE *fileOut3= fopen("greenW1.dat","w");
  if(fileOut3 == NULL) {printf("file %s not found\n", "greenW1.dat"); exit(1);}
  writeToFile_cMatrixFunction(fileOut3, &g0_matsubara, &model);
  fclose(fileOut3);  
  
  
  free_dMatrixFunction(&g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  free_cMatrixFunction(&g0_matsubara_Read);
  free_Model(&model);
  
  
  return Nerror;
}


FILE * fopenSafe(char fileName[], char mode[]){
  FILE * file = fopen(fileName, mode);
  if(file == NULL) {
    printf("error: file %s not found.\nterminated.\n", fileName); 
    exit(1);
  }
  else {
    printf("\nopening file %s\n", fileName);
  }
  return file;
}


int test_dmft(int verbose) {
  int Nerror=0;
  FILE * fileModel = fopenSafe("testInputFiles/dmft.model","rt");
  Model model;
  read_Model(fileModel,&model);
  fclose(fileModel);
  
  cMatrixFunction g0_matsubara;
  dMatrixFunction g0_tau;
  init_cMatrixFunction(&g0_matsubara, &model);
  init_dMatrixFunction(&g0_tau, &model);
  
  FILE * fileModel = fopenSafe("","rt");
  Model model;
  read_Model(fileModel,&model);
  fclose(fileModel);
  
  
  calculate_G0_matsubara(&g0_matsubara, &model, NULL);
  FILE * fileOut2 = fopenSafe("green_iwn.dat","w");
  writeToFile_cMatrixFunction(fileOut2, &g0_matsubara, &model);
  fclose(fileOut2);  

  calculate_G0_tau(&g0_matsubara,&g0_tau);
  
  FILE * fileOut = fopenSafe("green_tau.dat","w");
  writeToFile_dMatrixFunction(fileOut, &g0_tau, &model);
  fclose(fileOut);  
  
  
  
  free_dMatrixFunction(&g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  free_Model(&model);
  
  return Nerror;
}


