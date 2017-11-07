#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../findGreenSymmetries/findGreenSymmetries.h"
#include "../oneBodyMatrix/oneBodyMatrix.h"
//#include "../util/stringUtil.h"
//#include "../util/utilities.h"
//#include "../util/arrays/array.h"



typedef struct {
  HoppingMatrix tMat;
  MultiplePositions sites;
  MultiplePositions superlattice;
  dMatrix hybFM;
  Symmetries sym; 
  GreenSymmetriesMatrix greenSymMat;
  unsigned int nSites;
  double beta;
  double mu;
  double muAux;
  double U;
  double auxU;
  double delta;
  int nUpdates;
  int nThermUpdates; // Therm for thermalization
  int cleanUpdate_i;
  int measure_i;
} Model;


void read_Model(FILE * fileModel, FILE * fileParams, Model * model) {

  init_HoppingMatrix(&model->tMat);
  init_MultiplePositions(&model->sites);
  init_MultiplePositions(&model->superlattice);
     
  readOperators_HoppingMatrix(fileModel, &model->tMat);
  readSites(fileModel, &model->sites, "sites");
  print_MultiplePositions(&model->sites,"sites");
  model->nSites = model->sites.n;
  
  readSites(fileModel, &model->superlattice, "superlattice");
  print_MultiplePositions(&model->superlattice,"superlattice");
  assert(model->superlattice.n == 3);
  
  defineSparse_HoppingMatrix(&model->tMat, &model->sites, &model->superlattice);
  print_HoppingMatrix(&model->tMat);
  
  init_dMatrix(&model->hybFM,model->nSites);
  calculate_hybFirstMoments(&model->tMat, &model->hybFM);
  
  
  int nSym=countLineFlag(fileModel, "symmetry_generators");
  printf("nsym=%d\n", nSym);
  initSymmetries(&model->sym, nSym, model->nSites);
  readSymmetries(fileModel, model->nSites, &model->sym, "symmetry_generators");
  printSymmetries(&model->sym);
  
  initGreenSymmetriesMatrix(&model->greenSymMat,model->nSites,&model->sym);
  printf("\nsymmetrized matrix:\n");
  printGreenSymmetriesMatrix(&model->greenSymMat);

  
  printf("\n");
  
  readDouble(fileParams, "U",     &model->U);
  readDouble(fileParams, "mu",    &model->mu);
  readDouble(fileParams, "beta",  &model->beta);  
  readDouble(fileParams, "delta", &model->delta);  
  
  readInt(fileParams, "nUpdates",      &model->nUpdates);  
  readInt(fileParams, "nThermUpdates", &model->nThermUpdates);  
  readInt(fileParams, "cleanUpdate_i", &model->cleanUpdate_i);  
  readInt(fileParams, "measure_i",     &model->measure_i);  
  
  model->muAux = model->mu - model->U/2.0;
  model->auxU = model->U/2.0;
  
}


void free_Model(Model * model) {
  freeGreenSymmetriesMatrix(&model->greenSymMat);
  freeSymmetries(&model->sym);
  
  free_dMatrix(&model->hybFM);
  free_HoppingMatrix(&model->tMat);
  free_MultiplePositions(&model->sites);
  free_MultiplePositions(&model->superlattice);
}





