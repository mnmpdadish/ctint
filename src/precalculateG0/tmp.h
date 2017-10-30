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

#define N_PTS 1000

typedef struct {
  HoppingMatrix tMat;
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

  init_HoppingMatrix(&model->tMat);
  init_MultiplePositions(&model->sites);
  init_MultiplePositions(&model->superlattice);
     
  readOperators_HoppingMatrix(file, &model->tMat);
  readSites(file, &model->sites, "sites");
  print_MultiplePositions(&model->sites,"sites");
  
  readSites(file, &model->superlattice, "superlattice");
  print_MultiplePositions(&model->superlattice,"superlattice");
  assert(model->superlattice.n == 3);
  
  defineSparse_HoppingMatrix(&model->tMat, &model->sites, &model->superlattice);
  print_HoppingMatrix(&model->tMat);
  
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
  free_HoppingMatrix(&model->tMat);
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
  cMatrix **mat;
  cMatrix *M0; //zeroth moment
  cMatrix *M1; //first moment
  cMatrix *M2; //second moment
  cMatrix *M3; //third moment
  double beta;
} cMatrixFunction;

void init_cMatrixFunction(cMatrixFunction * cMatFun, Model * model) {

  printf("salut\n"); fflush(stdout);
  cMatFun->M0 = (cMatrix *) malloc(N_PTS * sizeof (cMatrix ));
  cMatFun->M1 = (cMatrix *) malloc(N_PTS * sizeof (cMatrix ));
  cMatFun->M2 = (cMatrix *) malloc(N_PTS * sizeof (cMatrix ));
  cMatFun->M3 = (cMatrix *) malloc(N_PTS * sizeof (cMatrix ));
  
  init_cMatrix(cMatFun->M0,model->sites.n);
  init_cMatrix(cMatFun->M1,model->sites.n);
  init_cMatrix(cMatFun->M2,model->sites.n);
  init_cMatrix(cMatFun->M3,model->sites.n);
  
  printf("salut\n"); fflush(stdout);
  cMatFun->mat = (cMatrix **) malloc(N_PTS * sizeof (cMatrix *));
  printf("salut\n"); fflush(stdout);
  
  int n;
  for(n=0; n<N_PTS; n++){
    cMatFun->mat[n] = (cMatrix *) malloc(N_PTS * sizeof (cMatrix));
    printf("%d ",n); fflush(stdout);
    init_cMatrix(cMatFun->mat[n],model->sites.n);
    reset_cMatrix(cMatFun->mat[n]);
  }
  
  /*
  //naming the n functions:
  int k;
  for(k=0; k<indepFun->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    nameGreenSymmetriesElement(&model->greenSymMat, i, j, indepFun->functions[k].name);
    indepFun->M0[k] = 0.0;
    indepFun->M1[k] = 0.0;
    indepFun->M2[k] = 0.0;
    indepFun->M3[k] = 0.0;
  }*/
  cMatFun->beta = model->beta;
}

void free_cMatrixFunction(cMatrixFunction * cMatFun) {
  free_cMatrix(cMatFun->M0);
  free_cMatrix(cMatFun->M1);
  free_cMatrix(cMatFun->M2);
  free_cMatrix(cMatFun->M3);
  int n;
  for(n=0; n<N_PTS; n++) free_cMatrix(cMatFun->mat[n]);
}

void calculateG0_matsubara(cMatrixFunction * g0_matsubara, Model * model) {
  int i,n;
  cMatrix tLoc;
  init_cMatrix(&tLoc,model->sites.n);
  calculate_HoppingMatrixLoc(&model->tMat, &tLoc);
  
  for(n=0; n<N_PTS; n++){
    reset_cMatrix(g0_matsubara->mat[n]);
    double complex z = I*(2.*n+1)*M_PI/g0_matsubara->beta;
    for(i=0; i<model->sites.n; i++) {
      printf("i=%d \n",i);
      printf("%f %f\n", creal(z + model->muAux), cimag(z + model->muAux));
      ELEM(g0_matsubara->mat[n], i, i) = z + model->muAux; // g = diagonal(muAux)
    }
    //cMatrixMatrixAddition(g0_matsubara->mat[n],&tLoc,g0_matsubara->mat[n], -1.0);  // g = g-tLoc
    //cMatrixMatrixAddition(&g0_matsubara->mat[n],&hyb_matsubara,&g0_matsubara->mat[n], 1.0); 
    print_cMatrix(g0_matsubara->mat[n]);
    print_cMatrix(&tLoc);
    invert_cMatrix(g0_matsubara->mat[n]); // g=g^-1
  }
}

/*
void calculateInversFourierTransform(cMatrixFunction * g0_matsubara, cMatrixFunction * g0_tau) {
  int k, n1, n2;
  for(k=0; k<g0_mat->n; k++){
    for(n1=0; n1<N_PTS; n1++){
      g0_tau->functions[k].data[n1] = 0.0;
      double tau = g0_mat->beta*n1/(N_PTS - 1);

      for(n2=N_PTS-1; n2>=0; n2--){
        double omega_n2 = (2.*n2+1)*M_PI/g0_mat->beta;
        double complex expFactor = 2.*cexp(-I*omega_n2*tau)/(g0_mat->beta) ;
        g0_tau->functions[k].data[n1] += expFactor * g0_mat->functions[k].data[n2];
      }
    }
  }
}
*/

/*
void calculateIndependant_G0_mat(IndepFunctionComplex * g0_mat, Model * model) {
  int k,n;
  for(k=0; k<g0_mat->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    double tLoc_ij = calculate_HoppingMatrixLoc_ij(&model->tMat, i, j);
    double a;
    if(i==j) a = -model->muAux + tLoc_ij; // only add mu on diagonal
    else a = tLoc_ij;
    g0_mat->M1[k] = 1.0;
    g0_mat->M2[k] = a;
    g0_mat->M3[k] = a*a; //+hybFM;
    for(n=0; n<N_PTS; n++){
      double complex z = I*(2.*n+1)*M_PI/g0_mat->beta;
      double complex hyb = 0.0;
      g0_mat->functions[k].data[n] = 1./( z - a - hyb);
      //if(n<100) printf("% 3.2e ",creal(g0->functions[k].data[n]));
    }
    printf("\n\n");
  }
}


void calculateInversFourierTransform(IndepFunctionComplex * g0_mat, IndepFunctionComplex * g0_tau) {
  int k, n1, n2;
  for(k=0; k<g0_mat->n; k++){
    for(n1=0; n1<N_PTS; n1++){
      g0_tau->functions[k].data[n1] = 0.0;
      double tau = g0_mat->beta*n1/(N_PTS - 1);

      for(n2=N_PTS-1; n2>=0; n2--){
        double omega_n2 = (2.*n2+1)*M_PI/g0_mat->beta;
        double complex expFactor = 2.*cexp(-I*omega_n2*tau)/(g0_mat->beta) ;
        g0_tau->functions[k].data[n1] += expFactor * g0_mat->functions[k].data[n2];
      }
    }
  }
}

void removeMoments_G0_mat(IndepFunctionComplex * g0_mat) { // can be in mat(tsubara) or (imaginary time) tau, as long a moments are defined
  int k, n;
  for(n=0; n<N_PTS; n++){
    double complex z = I*(2.*n+1)*M_PI/g0_mat->beta;
    for(k=0; k<g0_mat->n; k++){
      //g0_mat->functions[k].data[n] -= g0_mat->M0[k];
      g0_mat->functions[k].data[n] -= g0_mat->M1[k]/z;
      g0_mat->functions[k].data[n] -= g0_mat->M2[k]/(z*z);
      g0_mat->functions[k].data[n] -= g0_mat->M3[k]/(z*z*z);
    }
  }
}

void addMoments_G0_tau(IndepFunctionComplex * g0_mat, IndepFunctionComplex * g0_tau) { // can be in mat(tsubara) or (imaginary time) tau, as long a moments are defined
  int k, n;
  double beta = g0_tau->beta;
  for(n=0; n<N_PTS; n++){
    double tau = g0_tau->beta*n/(N_PTS - 1);
    for(k=0; k<g0_tau->n; k++){  
      //g0_tau->functions[k].data[n] =0.0;
      g0_tau->functions[k].data[n] += -g0_mat->M1[k]*0.5;
      g0_tau->functions[k].data[n] +=  g0_mat->M2[k]*0.5*(tau-0.5*beta);
      g0_tau->functions[k].data[n] +=  g0_mat->M3[k]*0.25*(beta-tau)*tau;
    }
  }
}


void calculateIndependant_G0_tau(IndepFunctionComplex * g0_mat, IndepFunctionComplex * g0_tau) {
  g0_tau->beta = g0_mat->beta;
  removeMoments_G0_mat(g0_mat);
  calculateInversFourierTransform(g0_mat, g0_tau);
  addMoments_G0_tau(g0_mat, g0_tau);
}
*/



void writeToFile_cMatrixFunction(FILE *fileOut, cMatrixFunction * cMatFun, Model * model) {
  int k,n;
  
  fprintf(fileOut, "# w_matsubara");
  for(k=0; k<model->greenSymMat.nIndep; k++) {
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    char nameReal[2];
    char nameImag[2];
    nameGreenSymmetriesElement(&model->greenSymMat, i, j, nameReal);
    nameGreenSymmetriesElement(&model->greenSymMat, i, j, nameImag);
    fprintf(fileOut, "          %s_re         %s_im", nameReal, nameImag);
  }
  
  for(n=0; n<N_PTS; n++){
    double omega_n = (2.*n+1)*M_PI/cMatFun->beta;
    fprintf(fileOut, "% 3.6e  ", omega_n);
    for(k=0; k<model->greenSymMat.nIndep; k++) {
      int i=model->greenSymMat.iFirstIndep[k];
      int j=model->greenSymMat.jFirstIndep[k];
      fprintf(fileOut, "% 3.6e % 3.6e  ", creal(ELEM(cMatFun->mat[n], i, j)), cimag(ELEM(cMatFun->mat[n], i, j)) );
    }
    fprintf(fileOut,"\n");
  }
}




