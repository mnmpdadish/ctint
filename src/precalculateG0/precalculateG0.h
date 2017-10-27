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
  double * M0; //zeroth moment
  double * M1; //first moment
  double * M2; //second moment
  double * M3; //third moment
  double beta;
} IndepFunctionComplex;

void init_IndepFunctionComplex(IndepFunctionComplex * indepFun, Model * model) {
  indepFun->n         = model->greenSymMat.nIndep;
  //indepFun->tau       = (functionComplex *) malloc(indepFun->n * sizeof (functionComplex));
  indepFun->functions = (functionComplex *) malloc(indepFun->n * sizeof (functionComplex));
  indepFun->M0 = (double *) malloc(indepFun->n * sizeof (double));
  indepFun->M1 = (double *) malloc(indepFun->n * sizeof (double));
  indepFun->M2 = (double *) malloc(indepFun->n * sizeof (double));
  indepFun->M3 = (double *) malloc(indepFun->n * sizeof (double));
  
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
  }
  indepFun->beta = model->beta;
}

void free_IndepFunctionComplex(IndepFunctionComplex * indepFun) {
  //free(indepFun->tau);
  free(indepFun->functions);
  free(indepFun->M0);
  free(indepFun->M1);
  free(indepFun->M2);
  free(indepFun->M3);
}

void calculateIndependant_G0_mat(IndepFunctionComplex * g0_mat, Model * model) {
  int k,n;
  for(k=0; k<g0_mat->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
    double tLoc_ij = calculate_tMatrixLoc_ij(&model->tMat, i, j);
    double a = -model->muAux + tLoc_ij;
    g0_mat->M1[k] = 1.0;
    g0_mat->M2[k] = a;
    g0_mat->M3[k] = a*a; //+hybFM;
    for(n=0; n<N_PTS; n++){
      double complex z = I*(2.*n+1)*M_PI/g0_mat->beta;
      double complex hyb = 0.0;
      g0_mat->functions[k].data[n] = 1./( z + model->muAux - tLoc_ij - hyb);
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

      for(n2=0; n2<N_PTS; n2++){
        double omega_n2 = (2.*n2+1)*M_PI/g0_mat->beta;
        double complex expFactor = (cos(-omega_n2*tau)+I*sin(-omega_n2*tau))/(g0_mat->beta) ;
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
      g0_mat->functions[k].data[n] -= g0_mat->M0[k];
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



/*
void IFT(const std::vector<_SiteVector>& bufferMAT, std::vector<_SiteVector>& bufferIT) {
		int tauN = 0;
		for(typename std::vector<_SiteVector>::iterator itIT = bufferIT.begin(); itIT != bufferIT.end(); itIT++, tauN++) {
			double tau = beta_*static_cast<double>(tauN)/static_cast<double>(bufferIT.size() - 1);
			itIT->clear();
			
			int frequencyMAT = 2*bufferMAT.size() - 1;
			for(typename std::vector<_SiteVector>::const_reverse_iterator itMAT = bufferMAT.rbegin(); itMAT != bufferMAT.rend(); itMAT++, frequencyMAT -= 2) {
				double omega = static_cast<double>(frequencyMAT)*M_PI/beta_;
				
				complex fact(2.*std::cos(-omega*tau)/beta_, 2.*std::sin(-omega*tau)/beta_);
				for(_Site i = 0; i < _SiteVector::VEC_DIM; ++i) 
					(*itIT)(i) += fact*(*itMAT)(i);
			}
		}
	};
	*/


void writeToFile_IndepFunctionComplex(FILE *fileOut, IndepFunctionComplex * indepFun) {
  int k,n;
  fprintf(fileOut, "# w_matsubara");
  for(k=0; k<indepFun->n; k++){
    fprintf(fileOut, "          %s_re         %s_im", (indepFun->functions[k].name ), (indepFun->functions[k].name ));
  }
  fprintf(fileOut, "\n");
      
  for(n=0; n<N_PTS; n++){
    double omega_n = (2.*n+1)*M_PI/indepFun->beta;
    //double complex z = I*omega_n;
    fprintf(fileOut, "% 3.6e  ", omega_n);
    //else break;
    for(k=0; k<indepFun->n; k++){
      fprintf(fileOut, "% 3.6e % 3.6e  ", creal(indepFun->functions[k].data[n]), cimag(indepFun->functions[k].data[n]));
      //moments: (they fit)
      //if(n<100) fprintf(fileOut, "% 3.6e % 3.6e  ", creal(indepFun->M2[k]/(z*z)), cimag(indepFun->M1[k]/z+indepFun->M3[k]/(z*z*z)) );
      //if(n<100) fprintf(fileOut, "% 3.6e % 3.6e  ", creal(indepFun->functions[k].data[n]-indepFun->M2[k]/(z*z)), cimag(indepFun->functions[k].data[n]-(indepFun->M1[k]/z+indepFun->M3[k]/(z*z*z))) );
    }
    fprintf(fileOut,"\n");
  }
}




