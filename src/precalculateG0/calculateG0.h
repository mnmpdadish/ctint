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

#define N_PTS_MAT 200
#define N_PTS_TAU 1000


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
  double auxU;
  double delta;
  int nUpdates;
  int cleanUpdate;
  int measureUpdate;
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
  readDouble(file, "delta",  &model->delta);  
  
  readInt(file, "nUpdates",  &model->nUpdates);  
  readInt(file, "cleanUpdate",  &model->cleanUpdate);  
  readInt(file, "measureUpdate",  &model->measureUpdate);  
  
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




// -------------------------------------------------------------------------

/*
typedef struct {
  double complex data[N_PTS];
  char name[2];
} functionComplex;
*/

typedef struct {
  cMatrix matrices[N_PTS_MAT];
  cMatrix M0; //zeroth moment
  cMatrix M1; //first moment
  cMatrix M2; //second moment
  cMatrix M3; //third moment
  double beta;
  unsigned int nAssigned;
} cMatrixFunction;

typedef struct {
  dMatrix matrices[N_PTS_TAU];
  double beta;
} dMatrixFunction;


void init_cMatrixFunction(cMatrixFunction * cMatFun, Model * model) {
  
  init_cMatrix(&cMatFun->M0,model->sites.n);
  init_cMatrix(&cMatFun->M1,model->sites.n);
  init_cMatrix(&cMatFun->M2,model->sites.n);
  init_cMatrix(&cMatFun->M3,model->sites.n);
  reset_cMatrix(&cMatFun->M0);
  reset_cMatrix(&cMatFun->M1);
  reset_cMatrix(&cMatFun->M2);
  reset_cMatrix(&cMatFun->M3);
  int n;
  for(n=0; n<N_PTS_MAT; n++){
    init_cMatrix(&cMatFun->matrices[n],model->sites.n);
    reset_cMatrix(&cMatFun->matrices[n]);
  }
  cMatFun->beta = model->beta;
}

void free_cMatrixFunction(cMatrixFunction * cMatFun) {
  free_cMatrix(&cMatFun->M0);
  free_cMatrix(&cMatFun->M1);
  free_cMatrix(&cMatFun->M2);
  free_cMatrix(&cMatFun->M3);
  int n;
  for(n=0; n<N_PTS_MAT; n++) free_cMatrix(&cMatFun->matrices[n]);
}


void init_dMatrixFunction(dMatrixFunction * dMatFun, Model * model) {
  int n;
  for(n=0; n<N_PTS_TAU; n++){
    init_dMatrix(&dMatFun->matrices[n],model->sites.n);
    reset_dMatrix(&dMatFun->matrices[n]);
  }
  dMatFun->beta = model->beta;
}

void free_dMatrixFunction(dMatrixFunction * dMatFun) {
  int n;
  for(n=0; n<N_PTS_TAU; n++) free_dMatrix(&dMatFun->matrices[n]);
}



void calculate_G0_matsubara(cMatrixFunction * g0_matsubara, Model * model, cMatrixFunction * hyb_matsubara) {
  int i,n;
  cMatrix tLoc;
  init_cMatrix(&tLoc,model->sites.n);
  calculate_HoppingMatrixLoc(&model->tMat, &tLoc);
  
  
  for(i=0; i<model->sites.n; i++) ELEM_VAL(g0_matsubara->M1, i, i) = 1.0; 
  printf("M1:\n"); print_cMatrix(&g0_matsubara->M1);
  for(i=0; i<model->sites.n; i++) ELEM_VAL(g0_matsubara->M2, i, i) = -model->muAux; 
  cMatrixMatrixAdditionInPlace(&g0_matsubara->M2, &tLoc, 1.0, 1.0);
  printf("M2:\n"); print_cMatrix(&g0_matsubara->M2);
  cMatrixMatrixMultiplication(&g0_matsubara->M2, &g0_matsubara->M2, &g0_matsubara->M3); //M2=M1*M1
  if(hyb_matsubara != NULL) 
    cMatrixMatrixAdditionInPlace(&g0_matsubara->M3,&hyb_matsubara->M1, 1.0, -1.0);  // g = g - hyb
  printf("M3:\n"); print_cMatrix(&g0_matsubara->M3);
  for(n=0; n<N_PTS_MAT; n++){
    reset_cMatrix(&g0_matsubara->matrices[n]);
    double complex z = I*(2.*n+1)*M_PI/g0_matsubara->beta;
    for(i=0; i<model->sites.n; i++) ELEM_VAL(g0_matsubara->matrices[n], i, i) = z + model->muAux; // g = diagonal(muAux)
    cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&tLoc, 1.0, -1.0);  // g = g - tLoc
    if(hyb_matsubara != NULL) 
      cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&hyb_matsubara->matrices[n], 1.0, -1.0);  // g = g - hyb
    invert_cMatrix(&g0_matsubara->matrices[n]); // g=g^-1
  }
}

//A=factorA*A+factorB*B
unsigned int dMatrix_cMatrixAdditionInPlace(dMatrix * A, cMatrix const * B, complex factorA, double complex factorB) {
  assert(A->N == B->N);
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=factorA*A->data[i]+creal(factorB*B->data[i]);
  return 0;
}

void calculateInversFourierTransform(cMatrixFunction * g0_matsubara, dMatrixFunction * g0_tau) {
  int n1, n2;
  for(n1=0; n1<N_PTS_TAU; n1++){
    //printf("salut3\n");
    reset_dMatrix(&g0_tau->matrices[n1]);
    double tau = g0_matsubara->beta*n1/(N_PTS_TAU - 1);
    for(n2=N_PTS_MAT-1; n2>=0; n2--){
      double omega_n2 = (2.*n2+1)*M_PI/g0_matsubara->beta;
      double complex expFactor = 2.*cexp(-I*omega_n2*tau)/(g0_matsubara->beta);
      //double complex expFactor2 = 1.*cexp(I*omega_n2*tau)/(g0_matsubara->beta);
      dMatrix_cMatrixAdditionInPlace(&g0_tau->matrices[n1],&g0_matsubara->matrices[n2], 1.0, expFactor);  // g = g -tLoc
      //dMatrix_cMatrixAdditionInPlace(&g0_tau->matrices[n1],&g0_matsubara->matrices[n2], 1.0, expFactor2);  // g = g -tLoc
      //printf("&g0_tau->matrices[n1] = \n");
      //print_cMatrix(&g0_tau->matrices[n1]);
    }
  }
}


void removeMoments_G0_matsubara(cMatrixFunction * g0_matsubara) {
  int n;
  for(n=0; n<N_PTS_MAT; n++){
    double complex z = I*(2.*n+1)*M_PI/g0_matsubara->beta;
    //cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&g0_matsubara->M0, 1.0, -1.0);  
    //printf("before:\n"); print_cMatrix(&g0_matsubara->matrices[n]);
    //printf("M1:\n"); print_cMatrix(&g0_matsubara->M1);
    cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&g0_matsubara->M1, 1.0, -1.0/z);  
    //printf("after: \n"); print_cMatrix(&g0_matsubara->matrices[n]);
    
    cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&g0_matsubara->M2, 1.0, -1.0/z/z);  
    cMatrixMatrixAdditionInPlace(&g0_matsubara->matrices[n],&g0_matsubara->M3, 1.0, -1.0/z/z/z);  
  }
}

void addMoments_G0_tau(cMatrixFunction * g0_matsubara, dMatrixFunction * g0_tau) { // can be in mat(tsubara) or (imaginary time) tau, as long a moments are defined
  int n;
  double beta = g0_tau->beta;
  for(n=0; n<N_PTS_TAU; n++){
    double tau = g0_tau->beta*n/(N_PTS_TAU - 1);
    dMatrix_cMatrixAdditionInPlace(&g0_tau->matrices[n],&g0_matsubara->M1, 1.0, -0.5);  
    dMatrix_cMatrixAdditionInPlace(&g0_tau->matrices[n],&g0_matsubara->M2, 1.0, 0.5*(tau-0.5*beta));  
    dMatrix_cMatrixAdditionInPlace(&g0_tau->matrices[n],&g0_matsubara->M3, 1.0, 0.25*(beta-tau)*tau);  
  }
}


void calculate_G0_tau(cMatrixFunction * g0_matsubara, dMatrixFunction * g0_tau) {
  g0_tau->beta = g0_matsubara->beta;
  removeMoments_G0_matsubara(g0_matsubara);
  calculateInversFourierTransform(g0_matsubara, g0_tau);
  addMoments_G0_tau(g0_matsubara, g0_tau);
}




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
  
  for(n=0; n<N_PTS_MAT; n++){
    fprintf(fileOut,"\n");
    double omega_n = (2.*n+1)*M_PI/cMatFun->beta;
    fprintf(fileOut, "% 3.6e  ", omega_n);
    for(k=0; k<model->greenSymMat.nIndep; k++) {
      int i=model->greenSymMat.iFirstIndep[k];
      int j=model->greenSymMat.jFirstIndep[k];
      fprintf(fileOut, "% 3.6e % 3.6e  ", creal(ELEM_VAL(cMatFun->matrices[n], i, j)), cimag(ELEM_VAL(cMatFun->matrices[n], i, j)) );
    }
  }
}



void writeToFile_dMatrixFunction(FILE *fileOut, dMatrixFunction * dMatFun, Model * model) {
  int k,n;
  
  fprintf(fileOut, "# tau        ");
  for(k=0; k<model->greenSymMat.nIndep; k++) {
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    char nameReal[2];
    nameGreenSymmetriesElement(&model->greenSymMat, i, j, nameReal);
    fprintf(fileOut, "          %s_re", nameReal);
  }
  
  for(n=0; n<N_PTS_TAU; n++){
    fprintf(fileOut,"\n");
    double tau = dMatFun->beta*n/(N_PTS_TAU - 1);
    fprintf(fileOut, "% 3.6e  ", tau);
    for(k=0; k<model->greenSymMat.nIndep; k++) {
      int i=model->greenSymMat.iFirstIndep[k];
      int j=model->greenSymMat.jFirstIndep[k];
      fprintf(fileOut, "% 3.6e  ", ELEM_VAL(dMatFun->matrices[n], i, j) );
    }
  }
}





void readFile_cMatrixFunction(FILE *fileIn, cMatrixFunction * cMatFun, Model * model) {
  rewind(fileIn);
  char tempbuff[2048];  //each line should not be above 2048 char long.
  int i=0, j=0, nRead;
  float w_matsubara;
  float valR, valI;
  char *pos;
  int n, n_matsubara=0;
  double complex indep_cMatrixValue[model->greenSymMat.nIndep];
  
        
  while(!feof(fileIn)) 
  {
    if (fgets(tempbuff,2048,fileIn)) {
      /*int lenString = strlen(tempbuff);
      if(tempbuff[lenString-1]!='\n' && tempbuff[lenString] != EOF) {
        printf("error. buffer too short? \\x%02x %d\n",tempbuff[lenString-1], lenString);
        printf("%s\n",tempbuff);
        printf("%d %d\n",tempbuff[lenString-1], 'a');
        
        exit(1);
      }*/
      unsigned int nElement = countElementInStr(tempbuff, " \t\n");
      if(tempbuff[0] == '#') continue;
      //else if(strBeginWithToken(tempbuff,"time")) sym->timeReversal=1;
      else if((tempbuff[0] != '\n') && nElement > 1) {
        pos = &tempbuff[0];
        valI=0.; valR=0.; w_matsubara=0.;
        nRead = sscanf(pos,"%e%n", &w_matsubara, &n);
        pos+=n;
        if(nRead!=1) {
          printf("Cannot read correctly the line: \n%s", tempbuff); 
          exit(1);
        }
        if(nElement > 2*model->greenSymMat.nIndep+1) {
          printf("Error, too many elements to read for the model\n");
          exit(1);
        }
        else if (nElement > 2*model->greenSymMat.nIndep+1) {
          printf("Error, not enough elements to read for the model\n");
          exit(1);
        }
        for(i=0; i<model->greenSymMat.nIndep; i++){
          nRead  = sscanf(pos,"%e%n", &valR, &n); pos+=n;
          nRead += sscanf(pos,"%e%n", &valI, &n); pos+=n;
          
          if(nRead!=2) {
            printf("Cannot read correctly the one-body line: \n%s", tempbuff); 
            exit(1);
          }
          indep_cMatrixValue[i] = valR + I*valI;
        }
        for(i=0;i<model->sites.n;i++) {
          for(j=0;j<model->sites.n;j++) {
            unsigned int index = model->greenSymMat.indexIndep[model->sites.n*i+j];
            ELEM_VAL(cMatFun->matrices[n_matsubara],i,j) = indep_cMatrixValue[index]; 
          }
        }
        //print_cMatrix(&cMatFun->matrices[n_matsubara]);
        n_matsubara++;
        if(n_mastubara >= N_PTS_MAT) break;
      }
      else break;
    }
    //printf("\n");
  }
  cMatFun->nAssigned=n_matsubara;
  return;
}

  
