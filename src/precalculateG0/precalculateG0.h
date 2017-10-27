#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../findGreenSymmetries/findGreenSymmetries.h"
#include "../oneBodyMatrix/oneBodyMatrix.h"
#include "../util/stringUtil.h"
#include "../util/utilities.h"
#include "../util/arrays/array.h"

#define N_PTS 2000

typedef struct {
  tMatrix tMat;
  MultiplePositions sites;
  MultiplePositions superlattice;
  dMatrix hybFM;
  Symmetries sym; 
  GreenSymmetriesMatrix greenSymMat;
  double beta;
  double mu;
  double muAux;
  double U;
} Model;


void read_Model(FILE * file, Model * model) {

  init_tMatrix(&model->tMat);
  init_MultiplePositions(&model->sites);
  init_MultiplePositions(&model->superlattice);
     
  readOperators_tMatrix(file, &model->tMat);
  readSites(file, &model->sites, "sites");
  print_MultiplePositions(&model->sites,"sites");
  
  readSites(file, &model->superlattice, "superlattice");
  print_MultiplePositions(&model->superlattice,"superlattice");
  assert(model->superlattice.n == 3);
  
  defineSparse_tMatrix(&model->tMat, &model->sites, &model->superlattice);
  print_tMatrix(&model->tMat);
  
  init_dMatrix(&model->hybFM,model->sites.n);
  calculate_hybFirstMoments(&model->tMat, &model->hybFM);
  
  
  int nSym=countLineFlag(file, "symmetry_generators");
  printf("nsym=%d\n", nSym);
  initSymmetries(&model->sym, nSym, model->sites.n);
  readSymmetries(file, model->sites.n, &model->sym, "symmetry_generators");
  printSymmetries(&model->sym);
  
  initGreenSymmetriesMatrix(&model->greenSymMat,model->sites.n,&model->sym);
  printf("\nsymmetrized matrix:\n");
  printGreenSymmetriesMatrix(&model->greenSymMat);
  
  printf("\n");
  
  
  readDouble(file, "U",     &model->U);
  readDouble(file, "mu",    &model->mu);
  readDouble(file, "beta",  &model->beta);  
  
  model->muAux = model->mu - model->U/2.0;
  
}


void free_Model(Model * model) {
  freeGreenSymmetriesMatrix(&model->greenSymMat);
  freeSymmetries(&model->sym);
  
  free_dMatrix(&model->hybFM);
  free_tMatrix(&model->tMat);
  free_MultiplePositions(&model->sites);
  free_MultiplePositions(&model->superlattice);
}




// -------------------------------------------------------------------------

typedef struct {
  double complex data[N_PTS];
  char name[2];
} functionComplex;


typedef struct {
  unsigned int n;
  //functionComplex * tau;
  functionComplex * functions;
  double * M1; //first moment
  double * M2; //second moment
  double * M3; //third moment
  double beta;
} IndepFunctionComplex;

void init_IndepFunctionComplex(IndepFunctionComplex * indepFun, Model * model) {
  indepFun->n         = model->greenSymMat.nIndep;
  //indepFun->tau       = (functionComplex *) malloc(indepFun->n * sizeof (functionComplex));
  indepFun->functions = (functionComplex *) malloc(indepFun->n * sizeof (functionComplex));
  indepFun->M1 = (double *) malloc(indepFun->n * sizeof (double));
  indepFun->M2 = (double *) malloc(indepFun->n * sizeof (double));
  indepFun->M3 = (double *) malloc(indepFun->n * sizeof (double));
  
  //naming the n functions:
  int k;
  for(k=0; k<indepFun->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    nameGreenSymmetriesElement(&model->greenSymMat, i, j, indepFun->functions[k].name);
    indepFun->M1[k] = 0.0;
    indepFun->M2[k] = 0.0;
    indepFun->M3[k] = 0.0;
  }
  indepFun->beta = model->beta;
}

void free_IndepFunctionComplex(IndepFunctionComplex * indepFun) {
  //free(indepFun->tau);
  free(indepFun->functions);
  free(indepFun->M1);
  free(indepFun->M2);
  free(indepFun->M3);
}

void precalculate_G0(IndepFunctionComplex * g0, Model * model) {
  int k,n;
  for(k=0; k<g0->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    double tLoc_ij = calculate_tMatrixLoc_ij(&model->tMat, i, j);
    double a = -model->muAux + tLoc_ij;
    g0->M1[k] = 1.0;
    g0->M2[k] = a;
    g0->M3[k] = a*a; //+hybFM;
    for(n=0; n<N_PTS; n++){
      double complex z = I*(2.*n+1)*M_PI/g0->beta;
      double complex hyb = 0.0;
      g0->functions[k].data[n] = 1./( z + model->muAux - tLoc_ij - hyb);
      //if(n<100) printf("% 3.2e ",creal(g0->functions[k].data[n]));
    }
    printf("\n\n");
  }
}


void writeToFile_IndepFunctionComplex(FILE *fileOut, IndepFunctionComplex * indepFun) {
  int k,n;
  fprintf(fileOut, "# w_matsubara");
  for(k=0; k<indepFun->n; k++){
    fprintf(fileOut, "          %s_re         %s_im", (indepFun->functions[k].name ), (indepFun->functions[k].name ));
  }
  fprintf(fileOut, "\n");
      
  for(n=0; n<N_PTS; n++){
    double omega_n = (2.*n+1)*M_PI/indepFun->beta;
    double complex z = I*omega_n;
    if(n<100) fprintf(fileOut, "% 3.6e  ", omega_n);
    else break;
    for(k=0; k<indepFun->n; k++){
      if(n<100) fprintf(fileOut, "% 3.6e % 3.6e  ", creal(indepFun->functions[k].data[n]), cimag(indepFun->functions[k].data[n]));
      if(n<100) fprintf(fileOut, "% 3.6e % 3.6e  ", creal(indepFun->M2[k]/(z*z)), cimag(indepFun->M1[k]/z+indepFun->M3[k]/(z*z*z)) );
    }
    fprintf(fileOut,"\n");
  }
}




