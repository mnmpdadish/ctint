#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "calculateG0/calculateG0.h"

#define VERTICES_BASIC_CAPACITY 256
#define N_BIN_TAU 400 //must be bigger than N_PTS_MAT: ideally 2*N_PTS_MAT
#define N_PTS_K 64


typedef struct {
  double tau;
  unsigned int site;
  int auxSpin;  // 1=up, 0=down.
} Vertex;

typedef struct {
  unsigned int N;
  unsigned int capacity;
  Vertex * m_vertex;
} VerticesChain;


typedef struct {
  double bin[N_BIN_TAU];
} dFunction;

typedef struct {
  unsigned int bin[N_BIN_TAU];
} iFunction;

typedef struct {
  double * indep_G_tau_sampled; //[mc->model.greenSymMat.nIndep]
  dMatrix G_tau0_sampled; 
  //dVector G_tau0_sampled; 
  //dVector G_tau0_sampled; 
  dFunction * indep_M_tau0_sampled;
  dFunction * indep_M_tau1_sampled;
  dFunction * indep_M_tau2_sampled;
  dFunction * indep_M_tau3_sampled;
  iFunction * nSamples;
  unsigned int nIndep;
} GreenAccumulator;


typedef struct {
  VerticesChain vertices;
  cMatrixFunction g0_matsubara; //non-intearcting Green Function.
  cMatrixFunction hyb_matsubara; //hybridisation
  
  cMatrix dummy1; 
  cMatrix dummy2; 
  
  dMatrixFunction g0_tau; //non-intearcting Green Function in imaginary time: to interpolate (in order to evaluate fast often).
  dMatrix *M_up;   // These two matrices are THE matrices
  dMatrix *M_down; // Same notation as Gull. 
  dMatrix *Mdummy_up; // We use pointer for the M matrices because we want to 
  dMatrix *Mdummy_down; // swap them at some point (by just swapping the pointer). Tested: it was 30% faster this way
  dVector Q, R;
  dVector Qtilde_up, Rtilde_up; //vectors to add
  dVector Qtilde_down, Rtilde_down; 
  double S_up,   Stilde_up; //corner to add
  double S_down, Stilde_down; 
  unsigned int nSites;
  Model model;
  double sign;
  //measurement values:
  double accumulated_sign;
  double accumulated_expOrder;
  cMatrixFunction accumulated_g_matsubara; //interacting Green Function to sample.
  dMatrix g_tau_0; //interacting Green Function at tau=0.
  //double complex expI[N_EXP];
  double KDirac; // same notation as Patrick, not sure about the math yet
  double density;
  GreenAccumulator g_tau_accumulator;
  
  unsigned long int nInsert;
  unsigned long int nRemove;
  unsigned long int nFlip;
} MonteCarlo;



void init_MonteCarlo(FILE * fileHyb, MonteCarlo * mc, Model * model) {
  
  init_dMatrixFunction(&mc->g0_tau, model);
  init_cMatrixFunction(&mc->accumulated_g_matsubara, model);
  init_cMatrixFunction(&mc->g0_matsubara, model);
  init_cMatrixFunction(&mc->hyb_matsubara, model);

  if(fileHyb!=NULL){  
    readFile_cMatrixFunction(fileHyb, &mc->hyb_matsubara, model);
    patch_HYB_matsubara(model, &mc->hyb_matsubara); //if nothing loaded, will patch with moments.
    calculate_G0_matsubara(&mc->g0_matsubara, model, &mc->hyb_matsubara, model->muAux, 0);
    
    cMatrixFunction g0_matsubara_tmp;
    init_cMatrixFunction(&g0_matsubara_tmp, model);
    calculate_G0_matsubara(&g0_matsubara_tmp, model, &mc->hyb_matsubara, model->muAux, 0);
    calculate_G0_tau(&g0_matsubara_tmp,&mc->g0_tau);
    free_cMatrixFunction(&g0_matsubara_tmp);
    
    //FILE * fileOut = fopenSafe("green_tau.dat","w",1);
    //writeToFile_dMatrixFunction(fileOut, &mc->g0_tau, model);
    //fclose(fileOut);  
  }
  else{
    calculate_G0_matsubara(&mc->g0_matsubara, model, NULL, model->muAux, 0);
  }
  
  //FILE * fileOut2 = fopenSafe("green_iwn.dat","w",1);
  //writeToFile_cMatrixFunction(fileOut2, &mc->g0_matsubara, model);
  //fclose(fileOut2);  

  
  init_cMatrix(&mc->dummy1, model->nSites);
  init_cMatrix(&mc->dummy2, model->nSites);
  init_dMatrix(&mc->g_tau_0, model->nSites);

  mc->M_up     = (dMatrix*)    malloc(sizeof(dMatrix));
  mc->M_down    = (dMatrix*)   malloc(sizeof(dMatrix));
  mc->Mdummy_up  = (dMatrix*)  malloc(sizeof(dMatrix));
  mc->Mdummy_down = (dMatrix*) malloc(sizeof(dMatrix));
  
  init_dMatrix(mc->M_up,0);
  init_dMatrix(mc->M_down,0);
  init_dMatrix(mc->Mdummy_up,0);
  init_dMatrix(mc->Mdummy_down,0);
  init_dVector(&mc->Q,0);
  init_dVector(&mc->R,0);
  init_dVector(&mc->Qtilde_up,0);
  init_dVector(&mc->Rtilde_up,0);
  init_dVector(&mc->Qtilde_down,0);
  init_dVector(&mc->Rtilde_down,0);
  
  mc->nSites = model->nSites;
  
  mc->g_tau_accumulator.nIndep = model->greenSymMat.nIndep;
  mc->g_tau_accumulator.indep_M_tau0_sampled = (dFunction*) malloc(model->greenSymMat.nIndep * sizeof(dFunction));
  mc->g_tau_accumulator.indep_M_tau1_sampled = (dFunction*) malloc(model->greenSymMat.nIndep * sizeof(dFunction));
  mc->g_tau_accumulator.indep_M_tau2_sampled = (dFunction*) malloc(model->greenSymMat.nIndep * sizeof(dFunction));
  mc->g_tau_accumulator.indep_M_tau3_sampled = (dFunction*) malloc(model->greenSymMat.nIndep * sizeof(dFunction));
  mc->g_tau_accumulator.indep_G_tau_sampled  = (double*)    malloc(model->greenSymMat.nIndep * sizeof(double));
  mc->g_tau_accumulator.nSamples             = (iFunction*) malloc(model->greenSymMat.nIndep * sizeof(iFunction));
  
  unsigned int k, n;
  for(k=0;k<model->greenSymMat.nIndep; k++){
    for(n=0;n<N_BIN_TAU; n++){
      mc->g_tau_accumulator.indep_M_tau0_sampled[k].bin[n] = 0.;
      mc->g_tau_accumulator.indep_M_tau1_sampled[k].bin[n] = 0.;
      mc->g_tau_accumulator.indep_M_tau2_sampled[k].bin[n] = 0.;
      mc->g_tau_accumulator.indep_M_tau3_sampled[k].bin[n] = 0.;
      mc->g_tau_accumulator.nSamples[k].bin[n] = 0;
    }
  }
  
  mc->vertices.capacity = VERTICES_BASIC_CAPACITY;
  mc->vertices.m_vertex = (Vertex*) malloc(mc->vertices.capacity * sizeof(Vertex));
  mc->vertices.N=0;
  mc->model=*model;
  mc->sign=1.0;
  mc->accumulated_sign=0;
  mc->accumulated_expOrder=0;
 
  mc->nRemove=0;
  mc->nInsert=0;
  mc->nFlip=0;
  //forcexp(-I*omega_n*tau);  
}

void free_MonteCarlo(MonteCarlo * mc) {
  free(mc->vertices.m_vertex);
  free_dMatrixFunction(&mc->g0_tau);
  free_cMatrixFunction(&mc->g0_matsubara);
  free_cMatrixFunction(&mc->accumulated_g_matsubara);

  free_dMatrix(mc->M_up);
  free_dMatrix(mc->M_down);
  free_dMatrix(mc->Mdummy_up);
  free_dMatrix(mc->Mdummy_down);
  free_dMatrix(&mc->Q);
  free_dMatrix(&mc->R);
  free_dMatrix(&mc->Qtilde_up);
  free_dMatrix(&mc->Rtilde_up);
  free_dMatrix(&mc->Qtilde_down);
  free_dMatrix(&mc->Rtilde_down);
  free(mc->g_tau_accumulator.indep_M_tau0_sampled);
  free(mc->g_tau_accumulator.indep_M_tau1_sampled);
  free(mc->g_tau_accumulator.indep_M_tau2_sampled);
  free(mc->g_tau_accumulator.indep_M_tau3_sampled);
  free(mc->g_tau_accumulator.indep_G_tau_sampled);
  free(mc->g_tau_accumulator.nSamples);
}

void Print_MonteCarlo(MonteCarlo * mc){
  int i, N=mc->vertices.N;
  printf("\n\nMonteCarlo state:\n");
  printf("\n"); for(i=0;i<N;i++) printf("%4.2f ",mc->vertices.m_vertex[i].tau); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].site); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].auxSpin); 
  /*printf("\nM_up:\n");
  print_dMatrix(mc->M_up);
  printf("\nM_down:\n");
  print_dMatrix(mc->M_down);
  */
  printf("\n");
  printf("\nM_up^-1:\n");
  copy_dMatrix(mc->M_up,mc->Mdummy_up);
  invert_dMatrix(mc->Mdummy_up);
  print_dMatrix(mc->Mdummy_up);
  printf("\nM_down^-1:\n");
  copy_dMatrix(mc->M_down,mc->Mdummy_down);
  invert_dMatrix(mc->Mdummy_down);
  print_dMatrix(mc->Mdummy_down);
  
}

double auxDown(Vertex const * vertex, MonteCarlo *mc) { return vertex->auxSpin ? -mc->model.delta : 1. + mc->model.delta;}
double auxUp(Vertex const * vertex, MonteCarlo *mc)   { return vertex->auxSpin ? 1. + mc->model.delta : -mc->model.delta;}
double urng() {return (double)rand()/(double)(RAND_MAX);}
unsigned int irng(unsigned int N) {return rand()%N;}


double green0(Vertex const * vertexI, Vertex const * vertexJ, MonteCarlo *mc) { 
  double diff_tau = vertexI->tau - vertexJ->tau + 1e-12; 

  double aps = 1.;
  if (diff_tau > .0) {
    diff_tau -= mc->model.beta;
    aps = -1.;
  }
  double ntau = fabs(diff_tau) * (N_PTS_TAU -1) /mc->model.beta;
  unsigned int ntau0= (unsigned int) ntau;
  
  double val_n  = ELEM_VAL(mc->g0_tau.matrices[ntau0],   vertexI->site, vertexJ->site);
  double val_n1 = ELEM_VAL(mc->g0_tau.matrices[ntau0+1], vertexI->site, vertexJ->site);
  
  return aps*((1. - (ntau - ntau0))*val_n + (ntau - ntau0)*val_n1);
}

/*
void swap_dMatrix_ptr(dMatrix *pA, dMatrix *pB){
  dMatrix *tmp = pA;
  pA = pB;
  pB = tmp;
}*///wut? does not work, keep the old way.

// add vertex at the end of the verticies:
int InsertVertex(MonteCarlo * mc) { 
  //printf("before insert\n"); fflush(stdout);
  Vertex newVertex;
  newVertex.tau = mc->model.beta*urng();
  newVertex.site= irng(mc->nSites);
  newVertex.auxSpin = (urng()<0.5)? 0 : 1;
  
  double S_up =   green0(&newVertex, &newVertex, mc) - auxUp(&newVertex, mc);
  double S_down = green0(&newVertex, &newVertex, mc) - auxDown(&newVertex, mc);
  //printf("S_up=%f\n",mc->S_up);
  //printf("S_down=%f\n",mc->S_down);
  
  if(mc->vertices.N>0){
    assert(mc->M_up->N == mc->vertices.N);
    assert(mc->M_down->N == mc->vertices.N);
    resize_dVector(&mc->R, mc->vertices.N);
    resize_dVector(&mc->Q, mc->vertices.N);
    resize_dVector(&mc->Rtilde_up, mc->vertices.N);
    resize_dVector(&mc->Rtilde_down, mc->vertices.N);
    resize_dVector(&mc->Qtilde_up, mc->vertices.N);
    resize_dVector(&mc->Qtilde_down, mc->vertices.N);
    
    reset_dVector(&mc->R);
    reset_dVector(&mc->Q);
    reset_dVector(&mc->Rtilde_up);
    reset_dVector(&mc->Rtilde_down);
    reset_dVector(&mc->Qtilde_up);
    reset_dVector(&mc->Qtilde_down);
    
    unsigned int n;
    for(n = 0; n < mc->vertices.N; n++) {
      Vertex vertexI = mc->vertices.m_vertex[n];
      mc->Q.data[n] = green0(&newVertex, &vertexI, mc);
      mc->R.data[n] = green0(&vertexI, &newVertex, mc);
    }
    
    dMatrixVectorProduct(mc->M_up,   &mc->Q, 1.0, &mc->Qtilde_up);
    dMatrixVectorProduct(mc->M_down, &mc->Q, 1.0, &mc->Qtilde_down);
    
    double Stilde_up =   1.0/( S_up   - dScalarProduct(&mc->R, &mc->Qtilde_up) );
    double Stilde_down = 1.0/( S_down - dScalarProduct(&mc->R, &mc->Qtilde_down) );
    
    double pAccept = -2.*mc->model.nSites*mc->model.beta*mc->model.auxU/ ( ((double)(mc->vertices.N+1))*Stilde_up*Stilde_down );
    //printf("pAcceptInsert=% 8.6f \n",pAccept);
    
    if(urng() < fabs(pAccept)) {
      //printf("before insert\n"); fflush(stdout);
      mc->nInsert++;
      //printf("% 9.5f ",pAccept); fflush(stdout);
			if(pAccept < .0) mc->sign *= -1;
			
//*  
      dVectorMatrixProduct(&mc->R, mc->M_up,   -Stilde_up, &mc->Rtilde_up);
      dVectorMatrixProduct(&mc->R, mc->M_down, -Stilde_down, &mc->Rtilde_down);
      
      scale_dVector(&mc->Qtilde_up,   -Stilde_up);
      scale_dVector(&mc->Qtilde_down, -Stilde_down);
      
      dAddRowColToInverse(mc->M_up,   &mc->Rtilde_up,   &mc->Qtilde_up,   Stilde_up,   mc->Mdummy_up);
      dAddRowColToInverse(mc->M_down, &mc->Rtilde_down, &mc->Qtilde_down, Stilde_down, mc->Mdummy_down);
      
      // swapping M with Mdummy:
      
      dMatrix *tmp = mc->M_up;
      mc->M_up = mc->Mdummy_up;
      mc->Mdummy_up = tmp;
      
      tmp = mc->M_down;
      mc->M_down = mc->Mdummy_down;
      mc->Mdummy_down = tmp;

      //printf("after insert\n"); fflush(stdout);
//*/
      
    }
    else return 0; //skip newVertex addition
  }
  else{
    double pAccept = -2.*mc->model.nSites*mc->model.beta*mc->model.auxU*S_up*S_down;
    //printf("%f\n",pAccept);
    if(urng() < fabs(pAccept)) {
      mc->nInsert++;
			if(pAccept < .0) mc->sign *= -1;
      resize_dMatrix(mc->M_up,1);
      resize_dMatrix(mc->M_down,1);
      
      mc->M_up->data[0]   = 1./S_up;
      mc->M_down->data[0] = 1./S_down;
    }
    else return 0; //skip newVertex addition
  }
  

  if(mc->vertices.N == mc->vertices.capacity) {
    mc->vertices.m_vertex = (Vertex*) realloc(mc->vertices.m_vertex, (mc->vertices.capacity *= 2) * sizeof(Vertex));
    printf("new size = %d \n",mc->vertices.capacity); fflush(stdout);
  }
  
  int N=mc->vertices.N;
  mc->vertices.m_vertex[N].tau  = newVertex.tau;
  mc->vertices.m_vertex[N].site = newVertex.site;
  mc->vertices.m_vertex[N].auxSpin = newVertex.auxSpin;
  
  mc->vertices.N++;
  
  //resize_dMatrix(mc->M_up, mc->vertices.N);
  //resize_dMatrix(mc->M_down, mc->vertices.N);
    
  
  //printf("after insert\n"); fflush(stdout);
  return 1;
}


// remove a vertex:
void RemoveVertex(MonteCarlo * mc) { 
  assert(mc->M_up->N == mc->vertices.N);
  assert(mc->M_down->N == mc->vertices.N);
  if(mc->vertices.N > 0){
    unsigned int p = irng(mc->vertices.N);
    double pAccept = (mc->vertices.N * ELEM(mc->M_up,p,p) * ELEM(mc->M_down,p,p)) / (-2.0*mc->model.nSites * mc->model.beta * mc->model.auxU);

    //printf("pAcceptRemove=% 8.6f \n",pAccept);
    if(urng() < fabs(pAccept)) {
      //printf("before remove\n"); fflush(stdout);
      mc->nRemove++;
			if(pAccept < .0) mc->sign *= -1;
      //printf("%d ==? %d\n",mc->M_up->N, mc->vertices.N);
      
      int N = mc->vertices.N-1;

      
      dMatrixSwapRows(mc->M_up,p,N);
      dMatrixSwapCols(mc->M_up,p,N);
      dMatrixSwapRows(mc->M_down,p,N);
      dMatrixSwapCols(mc->M_down,p,N);
      
      dSchurComplement(mc->M_up,  mc->Mdummy_up);
      dSchurComplement(mc->M_down,mc->Mdummy_down);
      
      // swapping M with Mdummy:
      dMatrix *tmp = mc->M_up;
      mc->M_up = mc->Mdummy_up;
      mc->Mdummy_up = tmp;
      
      tmp = mc->M_down;
      mc->M_down = mc->Mdummy_down;
      mc->Mdummy_down = tmp; 
      
      mc->vertices.N--;  // the vertex is not deleted, it is just forgotten
      
      mc->vertices.m_vertex[p].tau     = mc->vertices.m_vertex[N].tau;
      mc->vertices.m_vertex[p].site    = mc->vertices.m_vertex[N].site;
      mc->vertices.m_vertex[p].auxSpin = mc->vertices.m_vertex[N].auxSpin;
      //printf("after remove\n"); fflush(stdout);
      
      //resize_dMatrix(mc->M_up, mc->vertices.N);
      //resize_dMatrix(mc->M_down, mc->vertices.N);
  
      
    }

  }
}


// spin-flip a vertex:
void FlipVertex(MonteCarlo * mc) { 
  assert(mc->M_up->N == mc->vertices.N);
  assert(mc->M_down->N == mc->vertices.N);
  if(mc->vertices.N > 0){
    unsigned int p = irng(mc->vertices.N);
    
    Vertex vertexToFlip = mc->vertices.m_vertex[p];
    double flipAuxDown = auxDown(&vertexToFlip, mc) - auxUp(&vertexToFlip, mc);
    double flipAuxUp = auxUp(&vertexToFlip, mc) - auxDown(&vertexToFlip, mc);
    
    double factUp   = (1. + flipAuxUp * ELEM(mc->M_up,p,p));
    double factDown = (1. + flipAuxDown*ELEM(mc->M_down,p,p));
    
    double pAccept = factUp*factDown;
    //printf("pAccept spin-flip=%f\n", pAccept);
    if(urng() < fabs(pAccept)) {
      mc->nFlip++;
			if(pAccept < .0) mc->sign *= -1;
      
      double factorUp=-flipAuxUp/factUp;				
      dCopyRowIntoVector(mc->M_up, &mc->R, p); //automatic resizing
      dCopyColIntoVector(mc->M_up, &mc->Q, p); //automatic resizing
      dAddOneElementToInverse(mc->M_up, &mc->R, &mc->Q, factorUp);
      
      double factorDown=-flipAuxDown/factDown;				
      dCopyRowIntoVector(mc->M_down, &mc->R, p); //automatic resizing
      dCopyColIntoVector(mc->M_down, &mc->Q, p); //automatic resizing
      dAddOneElementToInverse(mc->M_down, &mc->R, &mc->Q, factorDown);
      
      mc->vertices.m_vertex[p].auxSpin = 1 - mc->vertices.m_vertex[p].auxSpin;  //if 0 become 1, if 1 become 0: flipped.
    }
  }
}


void CleanUpdate(MonteCarlo * mc) { 
  //printf("salut1\n"); fflush(stdout);
  resize_dMatrix(mc->Mdummy_up, mc->vertices.N);
  resize_dMatrix(mc->Mdummy_down, mc->vertices.N);
  //printf("salut2\n"); fflush(stdout);
  unsigned int i,j;
  for(i = 0; i < mc->vertices.N; i++) {
    Vertex vertexI = mc->vertices.m_vertex[i];
    for(j = 0; j < mc->vertices.N; j++) {
      Vertex vertexJ = mc->vertices.m_vertex[j];
      ELEM(mc->Mdummy_up,i,j) =   green0(&vertexJ, &vertexI, mc);
      ELEM(mc->Mdummy_down,i,j) = green0(&vertexJ, &vertexI, mc);
    }
    ELEM(mc->Mdummy_up,i,i) -= auxUp(&vertexI, mc);
    ELEM(mc->Mdummy_down,i,i) -= auxDown(&vertexI, mc);
  }
  //printf("salut3\n"); fflush(stdout);
  invert_dMatrix(mc->Mdummy_up);
  invert_dMatrix(mc->Mdummy_down);
  //printf("salut4\n"); fflush(stdout);
  
  // check what is the highest difference between Clean Update and current M matrices.
  double max_error=0.0;
  for(i = 0; i < mc->vertices.N; i++) {
    for(j = 0; j < mc->vertices.N; j++) {
      double error;
      error = fabs(ELEM(mc->Mdummy_up,i,j) - ELEM(mc->M_up,i,j));
      if(error > max_error ) max_error = error;
      error = fabs(ELEM(mc->Mdummy_down,i,j) - ELEM(mc->M_down,i,j));
      if(error > max_error ) max_error = error;
    }
  }
  printf("clean update max error= %4.5e\n",max_error);
  
  // swapping M with Mdummy:
      
  dMatrix *tmp = mc->M_up;
  mc->M_up = mc->Mdummy_up;
  mc->Mdummy_up = tmp;
  
  tmp = mc->M_down;
  mc->M_down = mc->Mdummy_down;
  mc->Mdummy_down = tmp;
      
}




int measure(MonteCarlo * mc) {
  unsigned int N = mc->vertices.N;
  unsigned int * sites;
  sites = (unsigned int*) malloc(N*sizeof(unsigned int));
    
  int p1,p2,k,i,j;
  double complex indepG_tau_sampled[mc->model.greenSymMat.nIndep];

  for(p1=0;p1<N;p1++) sites[p1] = mc->vertices.m_vertex[p1].site;

  const double DeltaInv = N_BIN_TAU/mc->model.beta;
  
  //printf("salut2\n");
  for(p1=0;p1<N;p1++) {
    for(p2=0;p2<N;p2++) {
      k = mc->model.greenSymMat.indexIndep[mc->model.nSites*sites[p1]+sites[p2]];
      double temp = mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2));
      double tau  = mc->vertices.m_vertex[p1].tau - mc->vertices.m_vertex[p2].tau; 
      if(tau < .0) {
        temp *= -1.;
        tau += mc->model.beta;
      }

      int index = DeltaInv*tau;
      double dTau = tau - ((double)index + 0.5)/DeltaInv;
      
      mc->g_tau_accumulator.indep_M_tau0_sampled[k].bin[index] += temp;
      temp*=dTau;
      mc->g_tau_accumulator.indep_M_tau1_sampled[k].bin[index] += temp;
      temp*=dTau;
      mc->g_tau_accumulator.indep_M_tau2_sampled[k].bin[index] += temp;
      temp*=dTau;
      mc->g_tau_accumulator.indep_M_tau3_sampled[k].bin[index] += temp;
      //mc->g_tau_accumulator.nSamples[k].bin[index]++;
    }
  }
  //printf("salut3\n"); fflush(stdout);
  
  //*
  for(k=0;k<mc->model.greenSymMat.nIndep;k++) {
    i = mc->model.greenSymMat.iFirstIndep[k];
    j = mc->model.greenSymMat.jFirstIndep[k];
    
    Vertex vertex0_I;
    vertex0_I.tau = 1e-13;
    vertex0_I.site= i;
    vertex0_I.auxSpin = 0;
    
    Vertex vertex0_J;
    vertex0_J.tau = 0;
    vertex0_J.site= j;
    vertex0_J.auxSpin = 0;
    
    indepG_tau_sampled[k] = mc->sign * green0(&vertex0_I, &vertex0_J, mc);
    for(p1=0;p1<N;p1++) {
      for(p2=0;p2<N;p2++) {
        //unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.nSites*sites[p1]+sites[p2]];
        //printf("salut=%d %d %d  %d %d\n", index, N*sites[p1]+sites[p2], mc->model.greenSymMat.nElement, sites[p1], sites[p2] ); fflush(stdout);
        indepG_tau_sampled[k] -= mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2)) * 
                                 green0(&vertex0_I,&mc->vertices.m_vertex[p2], mc) * 
                                 green0(&mc->vertices.m_vertex[p1],&vertex0_J, mc);
        
        if(p1==p2) {
          mc->KDirac += mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2));
          //printf("%f %f\n",green0(&vertex0_I,&mc->vertices.m_vertex[p1], mc),green0(&mc->vertices.m_vertex[p1],&vertex0_I, mc));
        }
      }
    }
    if(i==j) mc->density += indepG_tau_sampled[k]/mc->model.greenSymMat.numberOfSiteAssociated[k];
    mc->g_tau_accumulator.indep_G_tau_sampled[k] += indepG_tau_sampled[k];
  }
  //*/
  
  
  
  
  //printf("salut4\n"); fflush(stdout);
  
  mc->accumulated_sign +=mc->sign;
  mc->accumulated_expOrder +=N;
  free(sites);
  //printf("salut5\n"); fflush(stdout);
  
  return 1;
}


void calculate_G_matsubara_from_G_tau_accumulator(MonteCarlo *mc, cMatrixFunction *green_matsubara, unsigned int nSamples) {
  unsigned int n, i, j, k;
  double complex indep_M_matsubara_sampled[mc->model.greenSymMat.nIndep];
  const double dTau = mc->model.beta/N_BIN_TAU;
  for(n=0; n<N_PTS_MAT; n++) {
    double omega_n = M_PI*(2.*n + 1.)/mc->model.beta;
    double complex iomega_n = I*omega_n;
    double complex fact = cexp(iomega_n*dTau);
    double lambda = 2.*sin(omega_n*dTau/2.) / (dTau*omega_n*(1.-omega_n*omega_n*dTau*dTau/24.)*nSamples); //still does not understand this line.
    
    for(k=0;k<mc->model.greenSymMat.nIndep;k++){
      double complex temp_matsubara=0.;
      double complex exp_factor = cexp(iomega_n*dTau/2.)/(mc->model.greenSymMat.numberOfSiteAssociated[k]);  //watch out important factor!
      for(i=0;i<N_BIN_TAU;i++){
        double complex coeff = lambda*exp_factor;
        temp_matsubara += coeff * mc->g_tau_accumulator.indep_M_tau0_sampled[k].bin[i];
        temp_matsubara += coeff * mc->g_tau_accumulator.indep_M_tau1_sampled[k].bin[i] * iomega_n;
        temp_matsubara += coeff * mc->g_tau_accumulator.indep_M_tau2_sampled[k].bin[i] * iomega_n * iomega_n /2.;
        temp_matsubara += coeff * mc->g_tau_accumulator.indep_M_tau3_sampled[k].bin[i] * iomega_n * iomega_n * iomega_n /6.;
        
        exp_factor *= fact;
      }
      indep_M_matsubara_sampled[k] = temp_matsubara;
    }
    
    for(i=0;i<mc->model.nSites;i++) {
      for(j=0;j<mc->model.nSites;j++) {
        unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.nSites*i+j];
        ELEM_VAL(mc->dummy1,i,j) = indep_M_matsubara_sampled[index] / mc->model.beta; 
      }
    }
    
    // the three lines are: g = g0 - g0*dummy1*g0
    cMatrixMatrixMultiplication(&mc->g0_matsubara.matrices[n], &mc->dummy1, &mc->dummy2); // dummy2 = g0*dummy1
    cMatrixMatrixMultiplication(&mc->dummy2, &mc->g0_matsubara.matrices[n], &mc->dummy1); // dummy1 = dummy2*g0
    cMatrixMatrixAddition(&mc->dummy1, &mc->g0_matsubara.matrices[n], &green_matsubara->matrices[n], -1.0); // g = g0 - dummy1
  }
}


//self_or_hyb = diag(z+mu) - tLoc - green^-1 - hyb_or_self    (note that this equation is symmetric for self or hyb).
void extract_self_or_hyb_from_green(MonteCarlo *mc, cMatrixFunction *self_or_hyb, cMatrixFunction const *hyb_or_self, cMatrixFunction const *green) {
  unsigned int n, i;
  cMatrix tLoc;
  cMatrix inverted_green;
  init_cMatrix(&tLoc,mc->model.nSites);
  init_cMatrix(&inverted_green,mc->model.nSites);
  calculate_HoppingMatrixLoc(&mc->model.tMat, &tLoc);
  
  for(n=0;n<N_PTS_MAT;n++){
    double complex z = I*(2.*n+1)*M_PI/mc->model.beta; 
    reset_cMatrix(&self_or_hyb->matrices[n]);
    
    //copy_cMatrix(&green_matsubara.matrices[n],&inverted_green);
    copy_cMatrix(&green->matrices[n],&inverted_green);
    invert_cMatrix(&inverted_green);
    
    for(i=0; i<mc->model.nSites; i++) ELEM_VAL(self_or_hyb->matrices[n], i, i) = z + mc->model.mu;      // self = diag(z + mu) 
    cMatrixMatrixAdditionInPlace(&self_or_hyb->matrices[n], &tLoc, 1.0, -1.0);                          // self -= tLoc 
    cMatrixMatrixAdditionInPlace(&self_or_hyb->matrices[n], &inverted_green, 1.0, -1.0);                // self -= g^-1 
    //cMatrixMatrixAdditionInPlace(&self_matsubara.matrices[n], &mc->hyb_matsubara.matrices[n], 1.0, -1.0); // self -= hyb 
    cMatrixMatrixAdditionInPlace(&self_or_hyb->matrices[n], &hyb_or_self->matrices[n], 1.0, -1.0); // self -= hyb 
  }
  free_cMatrix(&tLoc);
  free_cMatrix(&inverted_green);
}

//g_lattice = sum_k  (1./(diag(z+mu) - t(k) - self)   
void integrate_green_lattice(MonteCarlo *mc, cMatrixFunction *green, cMatrixFunction const *self) {
  unsigned int n,i,j,k;
  cMatrix tK;
  cMatrix matrix_to_invert;
  init_cMatrix(&tK,mc->model.nSites);
  init_cMatrix(&matrix_to_invert,mc->model.nSites);

  for(n=0;n<N_PTS_MAT;n++) reset_cMatrix(&green->matrices[n]);
  
  //double factor = ((double)mc->nSites) / ( ((double)N_PTS_K)*((double)N_PTS_K) ); 
  double factor = 1.0 / ( ((double)N_PTS_K)*((double)N_PTS_K) ); 
  printf("factor=%f\n",factor);
  for(i=0;i<N_PTS_K;i++){
    //printf("%d\n",i);
    double kx = i*2*M_PI/N_PTS_K;
    for(j=0;j<N_PTS_K;j++){
      double ky = j*2*M_PI/N_PTS_K;
      calculate_tMatrixK_2D(&mc->model.tMat, &tK, kx, ky);
      for(n=0;n<N_PTS_MAT;n++){
        //printf("n=%d\n",n);
    
        double complex z = I*(2.*n+1)*M_PI/mc->model.beta; 
        reset_cMatrix(&matrix_to_invert);
        for(k=0; k<mc->model.nSites; k++) ELEM_VAL(matrix_to_invert, k, k) = z + mc->model.mu;      // g_i  = diag(z + mu) 
        cMatrixMatrixAdditionInPlace(&matrix_to_invert, &tK, 1.0, -1.0);                            // g_i -= t(k) 
        cMatrixMatrixAdditionInPlace(&matrix_to_invert, &self->matrices[n], 1.0, -1.0);             // g_i -= self 
        invert_cMatrix(&matrix_to_invert);                                                          // g_i  = g_i^-1 
        cMatrixMatrixAdditionInPlace(&green->matrices[n], &matrix_to_invert, 1.0, factor );         // g += g_i * Nc / N_k^2
      }
    }
  }
  free_cMatrix(&tK);
  free_cMatrix(&matrix_to_invert);
}



int outputMeasure(MonteCarlo * mc, unsigned int nSamples, unsigned long int iteration) {
  unsigned int n, i, j;
  
  cMatrixFunction self_matsubara;
  init_cMatrixFunction(&self_matsubara, &mc->model);
  
  double mu_new = mc->model.mu;
  
  if(iteration>0){
    
    for(i=0;i<mc->model.nSites;i++) {
      for(j=0;j<mc->model.nSites;j++) {
        unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.nSites*i+j];
        ELEM_VAL(mc->g_tau_0,i,j) += mc->g_tau_accumulator.indep_G_tau_sampled[index]/nSamples;
      }
    }
           
    double meas_sign = mc->accumulated_sign/nSamples;
    double meas_k = (mc->accumulated_expOrder/nSamples)/meas_sign;
    double meas_n = (2.0*mc->density/nSamples)/meas_sign;

    fappend("mu.dat",   "#    mu", iteration, mc->model.mu, 1);
    fappend("n.dat",    "#     n", iteration, meas_n, 1);
    fappend("sign.dat", "#  sign", iteration, meas_sign, 1);
    fappend("k.dat",    "#     k", iteration, meas_k, 1);
    
    double occupation =  -(mc->KDirac/(mc->model.beta*nSamples*mc->model.nSites) + (mc->model.muAux - mc->model.mu) )/mc->model.U;  
    printf("\nnSamples = %d\nsign = %f\nexpOrder =% f\nn =% f\nKDirac_n =% f\nnInsert = %lu\nnRemove = %lu\nnFlip = %lu\n",
           nSamples, meas_sign, meas_k, 2.0*mc->density/nSamples, 2.0*occupation,
           mc->nInsert, mc->nRemove, mc->nFlip);
    
    for(n=0; n<N_PTS_MAT; n++) scale_cMatrix(&mc->accumulated_g_matsubara.matrices[n],1.0/(nSamples*meas_sign));
  
    // extract Green:
    cMatrixFunction green_matsubara;
    init_cMatrixFunction(&green_matsubara, &mc->model);
  
    //------------------------------------------------
    calculate_G_matsubara_from_G_tau_accumulator(mc, &green_matsubara, nSamples);
    //------------------------------------------------
    char greenFileName[256];
    sprintf(greenFileName, "green%lu.dat", iteration); // puts string into buffer
    FILE *fileGreen = fopenSafe(greenFileName,  "w",1);
    writeToFile_cMatrixFunction(fileGreen, &green_matsubara, &mc->model);
    fclose(fileGreen);

    
    // extract self:
    //------------------------------------------------
    extract_self_or_hyb_from_green(mc, &self_matsubara, &mc->hyb_matsubara, &green_matsubara );
    //------------------------------------------------
    char selfFileName[256];
    sprintf(selfFileName, "self%lu.dat", iteration); // puts string into buffer
    FILE *fileSelf = fopenSafe(selfFileName,  "w",1);
    writeToFile_cMatrixFunction(fileSelf, &self_matsubara, &mc->model);
    fclose(fileSelf);
  }
  else if(iteration==0){ // for iteration 0, start from an empty self.
    for(n=0;n<N_PTS_MAT;n++) reset_cMatrix(&self_matsubara.matrices[n]);  // self = 0    
  }
  
  // integrate new green:
  cMatrixFunction new_green;
  init_cMatrixFunction(&new_green, &mc->model);
  
  //------------------------------------------------
  integrate_green_lattice(mc, &new_green, &self_matsubara);
  //------------------------------------------------
  
  // extract new hyb:
  cMatrixFunction new_hyb_matsubara;
  init_cMatrixFunction(&new_hyb_matsubara, &mc->model);
  //------------------------------------------------
  extract_self_or_hyb_from_green(mc, &new_hyb_matsubara, &self_matsubara, &new_green);
  //------------------------------------------------
  char hybFileName[256];
  sprintf(hybFileName, "hyb%lu.dat", iteration+1); // puts string into buffer
  FILE *fileHyb = fopenSafe(hybFileName, "w", 1);
  writeToFile_cMatrixFunction(fileHyb, &new_hyb_matsubara, &mc->model);
  fclose(fileHyb);
  
  
  char paramsN[256];
  char paramsNp1[256];
  char newMu[256];
  sprintf(paramsN,   "params%lu", iteration); // puts string into buffer
  sprintf(paramsNp1, "params%lu", iteration+1); // puts string into buffer
  sprintf(newMu, "%12.10f", mu_new); // puts string into buffer
  
  const char* paramName[] = { "mu" };
  const char* replacementValue[] = { newMu };
  copyFile(paramsN, paramsNp1, paramName, replacementValue, 1); 
  
  return 1;
}


