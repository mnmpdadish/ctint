#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "calculateG0.h"


int test_Plaquettes2x2(int verbose) {
  int Nerror=0;
  
  FILE * fileModel = fopenSafe("testInputFiles/plaquette2x2.model","rt",verbose);
  FILE * fileParams = fopenSafe("testInputFiles/params2x2","rt",verbose);
  Model model;
  read_Model(fileModel, fileParams, &model, verbose);
  fclose(fileModel);
  fclose(fileParams);
  
  cMatrixFunction g0_matsubara;
  cMatrixFunction g0_matsubara_Read;
  dMatrixFunction g0_tau;
  init_cMatrixFunction(&g0_matsubara, &model);
  init_cMatrixFunction(&g0_matsubara_Read, &model);
  init_dMatrixFunction(&g0_tau, &model);
  
  calculate_G0_matsubara(&g0_matsubara, &model, NULL, model.muAux, verbose);
  
  FILE *fileOut2= fopenSafe("greenW0.dat","w",verbose);
  writeToFile_cMatrixFunction(fileOut2, &g0_matsubara, &model);
  fclose(fileOut2);  


  calculate_G0_tau(&g0_matsubara,&g0_tau);
  
  FILE *fileOut = fopenSafe("green0.dat","w",verbose);
  writeToFile_dMatrixFunction(fileOut, &g0_tau, &model);
  fclose(fileOut);  
  
  printf("g0_tau.matrices[5].data[0]=%f\n", g0_tau.matrices[5].data[0]);
  printf("g0_tau.matrices[5].data[1]=%f\n", g0_tau.matrices[5].data[1]);
  printf("g0_tau.matrices[5].data[3]=%f\n", g0_tau.matrices[5].data[3]);
  if(!doubleEqual(g0_tau.matrices[5].data[0], -0.74174358)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[1],  0.13809594)) Nerror+=1;
  if(!doubleEqual(g0_tau.matrices[5].data[3],  0.19526082)) Nerror+=1;
  
  
  FILE *fileIn = fopenSafe("greenW0.dat","r",verbose);
  readFile_cMatrixFunction(fileIn, &g0_matsubara_Read, &model);
  fclose(fileIn);  
  
  FILE *fileOut3= fopenSafe("greenW1.dat","w",verbose);
  writeToFile_cMatrixFunction(fileOut3, &g0_matsubara, &model);
  fclose(fileOut3);  
  
  
  free_dMatrixFunction(&g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  free_cMatrixFunction(&g0_matsubara_Read);
  free_Model(&model);
  
  
  return Nerror;
}



int test_dmft(int verbose) {
  int Nerror=0;
  
  FILE * fileModel = fopenSafe("testInputFiles/dmft.model","rt",verbose);
  FILE * fileParams = fopenSafe("testInputFiles/paramsDmft","rt",verbose);
  Model model;
  read_Model(fileModel, fileParams, &model, verbose);
  fclose(fileModel);
  fclose(fileParams);
  
  cMatrixFunction g0_matsubara;
  dMatrixFunction g0_tau;
  cMatrixFunction hyb_matsubara;
  init_cMatrixFunction(&g0_matsubara, &model);
  init_dMatrixFunction(&g0_tau, &model);
  init_cMatrixFunction(&hyb_matsubara, &model);
  
  FILE * fileHyb = fopenSafe("testInputFiles/hyb1.dat","rt",verbose);
  readFile_cMatrixFunction(fileHyb, &hyb_matsubara, &model);
  patch_HYB_matsubara(&model, &hyb_matsubara);
  fclose(fileHyb);
  
  
  
  FILE * fileHybOut = fopenSafe("hybOut.dat","w",verbose);
  writeToFile_cMatrixFunction(fileHybOut, &hyb_matsubara, &model);
  fclose(fileHybOut);
  
  
  
  calculate_G0_matsubara(&g0_matsubara, &model, &hyb_matsubara, model.muAux, verbose);
  calculate_G0_tau(&g0_matsubara,&g0_tau);
  FILE * fileOut2 = fopenSafe("green_iwn.dat","w",verbose);
  writeToFile_cMatrixFunction(fileOut2, &g0_matsubara, &model);
  fclose(fileOut2);  

  
  FILE * fileOut = fopenSafe("green_tau.dat","w",verbose);
  writeToFile_dMatrixFunction(fileOut, &g0_tau, &model);
  fclose(fileOut);  
  
  free_dMatrixFunction(&g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  free_Model(&model);
  
  return Nerror;
}


