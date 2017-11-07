#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "calculateG0/calculateG0.h"

#define VERTICES_BASIC_CAPACITY 256
#define N_EXP 2000


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
  VerticesChain vertices;
  cMatrixFunction g0_matsubara; //non-intearcting Green Function.
  cMatrixFunction hyb_matsubara; //hybridisation
  
  cMatrix dummy1; 
  cMatrix dummy2; 
  
  dMatrixFunction g0_tau; //non-intearcting Green Function in imaginary time: to interpolate (in order to evaluate fast often).
  dMatrix *M_up;   // These two matrices are THE matrices
  dMatrix *M_down; // Same notation as Gull. 
  dMatrix *Mdummy_up; // We use pointer for the M matrices because we want to 
  dMatrix *Mdummy_down; // swap them at some point (by just swapping the pointer). Tested: it was 15% faster this way
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
  dMatrix accumulated_g_tau; //interacting Green Function to sample.
  double complex expI[N_EXP];
  double KDirac; // same notation as Patrick, not sure about the math yet
  double density;
} MonteCarlo;


void init_MonteCarlo(FILE * fileHyb, MonteCarlo * mc, Model * model) {
  init_dMatrixFunction(&mc->g0_tau, model);
  init_cMatrixFunction(&mc->accumulated_g_matsubara, model);
  init_cMatrixFunction(&mc->g0_matsubara, model);
  init_cMatrixFunction(&mc->hyb_matsubara, model);

  //printf("salut"); fflush(stdout);
  //FILE * fileHyb = fopenSafe("files/hyb99.dat","rt", 1);
  readFile_cMatrixFunction(fileHyb, &mc->hyb_matsubara, model);
  patch_HYB_matsubara(model, &mc->hyb_matsubara);
  fclose(fileHyb);
  
  calculate_G0_matsubara(&mc->g0_matsubara, model, &mc->hyb_matsubara, model->muAux, 1);
  
  init_cMatrix(&mc->dummy1, model->sites.n);
  init_cMatrix(&mc->dummy2, model->sites.n);
  init_dMatrix(&mc->accumulated_g_tau, model->sites.n);
  
  cMatrixFunction g0_matsubara_tmp;
  init_cMatrixFunction(&g0_matsubara_tmp, model);
  calculate_G0_matsubara(&g0_matsubara_tmp, model, &mc->hyb_matsubara, model->muAux, 0);
  calculate_G0_tau(&g0_matsubara_tmp,&mc->g0_tau);
  free_cMatrixFunction(&g0_matsubara_tmp);
  
  mc->M_up     = malloc(sizeof(dMatrix));
  mc->M_down    = malloc(sizeof(dMatrix));
  mc->Mdummy_up  = malloc(sizeof(dMatrix));
  mc->Mdummy_down = malloc(sizeof(dMatrix));
  
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
  
  mc->nSites = model->sites.n;
  
  mc->vertices.capacity = VERTICES_BASIC_CAPACITY;
  mc->vertices.m_vertex = malloc(mc->vertices.capacity * sizeof(Vertex));
  mc->vertices.N=0;
  mc->model=*model;
  mc->sign=1.0;
  mc->accumulated_sign=0;
  mc->accumulated_expOrder=0;
 
  unsigned int i;
  for(i=0;i<N_EXP;i++) mc->expI[i] = cexp(-I*2*M_PI*((double)i)/((double) N_EXP)); //precalculate every green function

  //forcexp(-I*omega_n*tau);  
}

void free_MonteCarlo(MonteCarlo * mc) {
  free(mc->vertices.m_vertex);
  free_dMatrixFunction(&mc->g0_tau);
  free_cMatrixFunction(&mc->g0_matsubara);
  free_cMatrixFunction(&mc->accumulated_g_matsubara);
  free_dMatrix(&mc->accumulated_g_tau);
  //free_dMatrixFunction(&mc->g0_tau_down);
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
}

void Print_MonteCarlo(MonteCarlo * mc){
  int i, N=mc->vertices.N;
  printf("MonteCarlo state:\n");
  printf("\n"); for(i=0;i<N;i++) printf("%4.2f ",mc->vertices.m_vertex[i].tau); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].site); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].auxSpin); 
  printf("\nM_up:\n");
  print_dMatrix(mc->M_up);
  printf("\nM_down:\n");
  print_dMatrix(mc->M_down);
  printf("\n");
}

double auxUp(Vertex const * vertex, MonteCarlo *mc)   { return vertex->auxSpin ? -mc->model.delta : 1. + mc->model.delta;}
double auxDown(Vertex const * vertex, MonteCarlo *mc) { return vertex->auxSpin ? 1. + mc->model.delta : -mc->model.delta;}
double urng() {return (double)rand()/(double)(RAND_MAX);}
unsigned int irng(unsigned int N) {return rand()%N;}


double green0(Vertex const * vertexI, Vertex const * vertexJ, MonteCarlo *mc) { 
  double diff_tau = vertexI->tau - vertexJ->tau +1e-13; 

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

	

// add vertex at the end of the verticies:
int InsertVertex(MonteCarlo * mc) { 
  
  Vertex newVertex;
  newVertex.tau = mc->model.beta*urng();
  newVertex.site= irng(mc->nSites);
  newVertex.auxSpin = irng(2);
  
  mc->S_up =   green0(&newVertex, &newVertex, mc) - auxUp(&newVertex, mc);
  mc->S_down = green0(&newVertex, &newVertex, mc) - auxDown(&newVertex, mc);
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
    unsigned int n;
    for(n = 0; n < mc->vertices.N; n++) {
      Vertex vertexI = mc->vertices.m_vertex[n];
      mc->R.data[n] = green0(&newVertex, &vertexI, mc);
      mc->Q.data[n] = green0(&vertexI, &newVertex, mc);
    }
			
    dMatrixVectorProduct(mc->M_up, &mc->Q, 1.0, &mc->Qtilde_up);
    dMatrixVectorProduct(mc->M_down, &mc->Q, 1.0, &mc->Qtilde_down);
    
    mc->Stilde_up = 1.0/( mc->S_up - dScalarProduct(&mc->R, &mc->Qtilde_up) );
    mc->Stilde_down = 1.0/( mc->S_down - dScalarProduct(&mc->R, &mc->Qtilde_down) );
    
    double pAcc = -2.*mc->model.sites.n*mc->model.beta*mc->model.auxU/ ( ((double)mc->vertices.N)*mc->Stilde_up*mc->Stilde_down );
    //printf("pAcc=%f\n",pAcc);
    if(urng() < fabs(pAcc)) {
			if(pAcc < .0) mc->sign *= -1;
      
      dVectorMatrixProduct(&mc->R, mc->M_up, -mc->Stilde_up, &mc->Rtilde_up);
      dVectorMatrixProduct(&mc->R, mc->M_down, -mc->Stilde_down, &mc->Rtilde_down);
      
      scale_dVector(&mc->Qtilde_up, -mc->Stilde_up);
      scale_dVector(&mc->Qtilde_down, -mc->Stilde_down);
      
      dAddRowColToInverse(mc->M_up,   &mc->Rtilde_up,   &mc->Qtilde_up,   mc->Stilde_up,   mc->Mdummy_up);
      dAddRowColToInverse(mc->M_down, &mc->Rtilde_down, &mc->Qtilde_down, mc->Stilde_down, mc->Mdummy_down);
      
      // swapping M with Mdummy:
      dMatrix *tmp = mc->M_up;
      mc->M_up = mc->Mdummy_up;
      mc->Mdummy_up = tmp;
      
      tmp = mc->M_down;
      mc->M_down = mc->Mdummy_down;
      mc->Mdummy_down = tmp;
      
      //old way: copy instead of swapping vector (slower):
      //copy_dMatrix(mc->Mdummy_up,mc->M_up);
      //copy_dMatrix(mc->Mdummy_down,mc->M_down);

      
    }
    else return 0; //skip newVertex addition
  }
  else{
    double pAcc = -2.*mc->model.sites.n*mc->model.beta*mc->model.auxU*mc->S_up*mc->S_down;
    //printf("%f\n",pAcc);
    if(urng() < fabs(pAcc)) {
      if(pAcc < .0) mc->sign *= -1;
      resize_dMatrix(mc->M_up,1);
      resize_dMatrix(mc->M_down,1);
      
      mc->M_up->data[0]   = 1./mc->S_up;
      mc->M_down->data[0] = 1./mc->S_down;
    }
    else return 0; //skip newVertex addition
  }
  
  if(mc->vertices.N == mc->vertices.capacity) {
    mc->vertices.m_vertex = realloc(mc->vertices.m_vertex, (mc->vertices.capacity *= 2) * sizeof(Vertex));
  }
  
  int N=mc->vertices.N;
  mc->vertices.m_vertex[N].tau  = newVertex.tau;
  mc->vertices.m_vertex[N].site = newVertex.site;
  mc->vertices.m_vertex[N].auxSpin = newVertex.auxSpin;
  
  mc->vertices.N++;
  
  return 1;
}




// remove a vertex:
void RemoveVertex(MonteCarlo * mc) { 
  assert(mc->M_up->N == mc->vertices.N);
  assert(mc->M_down->N == mc->vertices.N);
  if(mc->vertices.N > 0){
    unsigned int p = irng(mc->vertices.N);
    //double factUp = -1./ELEM(mc->M_up,p,p);
    //double factDown = -1./ELEM(mc->M_down,p,p);
		//double pAcc = mc->vertices.N / (-2.0*mc->model.sites.n * mc->model.beta * mc->model.auxU * factUp * factDown);
    double pAcc = (mc->vertices.N*ELEM(mc->M_up,p,p)*ELEM(mc->M_down,p,p) ) / (-2.0*mc->model.sites.n * mc->model.beta * mc->model.auxU);
    
    if(urng() < fabs(pAcc)) {
      if(pAcc < .0) mc->sign *= -1;
      //printf("%d ==? %d\n",mc->M_up->N, mc->vertices.N);
      
      mc->vertices.N--;  // the vertex is not deleted, it is just forgotten
      int N = mc->vertices.N;
      dMatrixSwapRows(mc->M_up,p,N);
      dMatrixSwapCols(mc->M_up,p,N);
      dMatrixSwapRows(mc->M_down,p,N);
      dMatrixSwapCols(mc->M_down,p,N);
      
      dSchurComplement(mc->M_up,  mc->Mdummy_up);
      dSchurComplement(mc->M_down,mc->Mdummy_down);
      
      copy_dMatrix(mc->Mdummy_up,  mc->M_up);
      copy_dMatrix(mc->Mdummy_down,mc->M_down);
      
      mc->vertices.m_vertex[p].tau     = mc->vertices.m_vertex[N].tau;
      mc->vertices.m_vertex[p].site    = mc->vertices.m_vertex[N].site;
      mc->vertices.m_vertex[p].auxSpin = mc->vertices.m_vertex[N].auxSpin;
    }
  }
}



void CleanUpdate(MonteCarlo * mc) { 
  assert(mc->M_up->N == mc->vertices.N);
  assert(mc->M_down->N == mc->vertices.N);
  unsigned int i,j;
  for(i = 0; i < mc->vertices.N; i++) {
    Vertex vertexI = mc->vertices.m_vertex[i];
    for(j = 0; j < mc->vertices.N; j++) {
      Vertex vertexJ = mc->vertices.m_vertex[j];
      ELEM(mc->M_up,i,j) = green0(&vertexI, &vertexJ, mc);
      ELEM(mc->M_down,i,j) = green0(&vertexI, &vertexJ, mc);
    }
    ELEM(mc->M_up,i,i) -= auxUp(&vertexI, mc);
    ELEM(mc->M_down,i,i) -= auxDown(&vertexI, mc);
  }
  invert_dMatrix(mc->M_up);
  invert_dMatrix(mc->M_down);
}












//#define N_PTS_OUT 200
// -------------------------------------------------------------------------

/*
typedef struct {
  double data[N_PTS_TAU];
} functionDouble;

typedef struct {
  unsigned int n;
  functionDouble * func;
} IndepFunctionDouble;

void init_IndepFunctionDouble(IndepFunctionDouble * indepFun, Model * model) {
  indepFun->n    = model->greenSymMat.nIndep;
  indepFun->func = (functionDouble *) malloc(indepFun->n * sizeof (functionDouble));

  int k;
  for(k=0; k<indepFun->n; k++){
    int i=model->greenSymMat.iFirstIndep[k];
    int j=model->greenSymMat.jFirstIndep[k];
  }
}

void free_IndepFunctionDouble(IndepFunctionDouble * indepFun) {
  free(indepFun->func);
}
*/


/*
double complex expI(double value, MonteCarlo * mc) { 
  //double diff_tau = vertexI->tau - vertexJ->tau; 

  //double aps = 1.;
  //if (value < .0) 
  double       nValue = fabs(diff_tau) * (N_EXP -1) /mc->model.beta;
  unsigned int ntau0= (unsigned int) ntau;
  
  double val_n  = ELEM_VAL(mc->g0_tau.matrices[ntau0], vertexI->site, vertexJ->site);
  double val_n1 = ELEM_VAL(mc->g0_tau.matrices[ntau0+1], vertexI->site, vertexJ->site);
  
  return aps*((1. - (ntau - ntau0))*val_n + (ntau - ntau0)*val_n1);
}
//


void build_dMatrix_from_indep(Model * model, dMatrix * matrix, double arrayIndep[]) {
  assert(matrix->N==model->sites.n);
  assert(matrix->N == sizeof(arrayIndep)/sizeof(arrayIndep[0]));
  unsigned int i,j, N=matrix->N;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      unsigned int index = model->greenSymMat.indexIndep[N*i+j];
      ELEM(matrix,i,j) += arrayIndep[index];
    }
  }
}

void build_cMatrix_from_indep(Model * model, cMatrix * matrix, double arrayIndep[]) {
  assert(matrix->N==model->sites.n);
  assert(matrix->N == sizeof(arrayIndep)/sizeof(arrayIndep[0]));
  unsigned int i,j, N=matrix->N;
  for(i=0;i<N;i++) {
    for(j=0;j<N;j++) {
      unsigned int index = model->greenSymMat.indexIndep[N*i+j];
      ELEM(matrix,i,j) += arrayIndep[index];
    }
  }
}

//*/

int measure(MonteCarlo * mc) {
  unsigned int N = mc->vertices.N;
  //double green0_tau1[N];
  //double green0_tau2[N];
  double complex exp_IomegaN_tau[N];
  double complex indepM_sampled[mc->model.greenSymMat.nIndep];
  double indepG_tau_sampled[mc->model.greenSymMat.nIndep];
  unsigned int sites[N];
  int n,p1,p2,k,i,j;
  
  
  for(k=0;k<mc->model.greenSymMat.nIndep;k++){
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
        //unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.sites.n*sites[p1]+sites[p2]];
        //printf("salut=%d %d %d  %d %d\n", index, N*sites[p1]+sites[p2], mc->model.greenSymMat.nElement, sites[p1], sites[p2] ); fflush(stdout);
        indepG_tau_sampled[k] -= mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2)) * 
                                 green0(&vertex0_I,&mc->vertices.m_vertex[p1], mc) * 
                                 green0(&mc->vertices.m_vertex[p2],&vertex0_J, mc);
        
        if(p1==p2) {
          mc->KDirac += mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2));
          //printf("%f %f\n",green0(&vertex0_I,&mc->vertices.m_vertex[p1], mc),green0(&mc->vertices.m_vertex[p1],&vertex0_I, mc));
        }
      }
    }
    if(i==j) mc->density += indepG_tau_sampled[k]*mc->model.greenSymMat.numberOfSiteAssociated[k] / mc->model.sites.n;
  }
  //*/
  
  for(i=0;i<mc->model.sites.n;i++) {
    for(j=0;j<mc->model.sites.n;j++) {
      unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.sites.n*i+j];
      ELEM_VAL(mc->accumulated_g_tau,i,j) += indepG_tau_sampled[index];
    }
  }
    
  
  //printf(" n= % 5.3f", indepG_tau_sampled[0]);
  /*for(p1=0;p1<N;p1++) {
    green0_tau1[p1] = green0(mc->vertices.m_vertex[p1]);
    exp_IomegaN_tau[p1] = cexp(-I*omega_n*tau); //precalculate every green function
  }//*/
    
  
  //double g_at_tau0;
  
  for(p1=0;p1<N;p1++) sites[p1] = mc->vertices.m_vertex[p1].site;
  
  //Print_MonteCarlo(mc);
  for(n=0;n<N_PTS_MAT;n++){
    double omega_n = (2.*n+1)*M_PI/mc->model.beta; 
    //reset_cMatrix(&mc->accumulated_g_matsubara.matrices[n]);
    reset_cMatrix(&mc->dummy1);
    reset_cMatrix(&mc->dummy2);
    for(p1=0;p1<N;p1++) {
      double tau = mc->vertices.m_vertex[p1].tau;
      exp_IomegaN_tau[p1] = cexp(-I*omega_n*tau); //precalculate every green function
    }
    for(k=0;k<mc->model.greenSymMat.nIndep; k++) indepM_sampled[k]=0.0;
    
    for(p1=0;p1<N;p1++) {
      for(p2=0;p2<N;p2++) {
        unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.sites.n*sites[p1]+sites[p2]];
        //printf("salut=%d %d %d  %d %d\n", index, N*sites[p1]+sites[p2], mc->model.greenSymMat.nElement, sites[p1], sites[p2] ); fflush(stdout);
        indepM_sampled[index] += exp_IomegaN_tau[p1]*conj(exp_IomegaN_tau[p2]) * mc->sign * 0.5 * ( ELEM(mc->M_up, p1, p2) + ELEM(mc->M_down, p1, p2));
      }
    }
    
    for(i=0;i<mc->model.sites.n;i++) {
      for(j=0;j<mc->model.sites.n;j++) {
        unsigned int index = mc->model.greenSymMat.indexIndep[mc->model.sites.n*i+j];
        ELEM_VAL(mc->dummy1,i,j) = indepM_sampled[index] / mc->model.beta; 
      }
    }
    
    // the three lines are: g = g0 - g0*dummy1*g0
    
    cMatrixMatrixMultiplication(&mc->g0_matsubara.matrices[n], &mc->dummy1, &mc->dummy2); // dummy2 = g0*dummy1
    cMatrixMatrixMultiplication(&mc->dummy2, &mc->g0_matsubara.matrices[n], &mc->dummy1); // dummy1 = dummy2*g0
    cMatrixMatrixAddition(&mc->dummy1, &mc->g0_matsubara.matrices[n], &mc->dummy2, -1.0); // dummy2 = g0 - dummy1
    
    /*
    if(n<3){
      printf("dummy %d ", n);
      print_cMatrix(&mc->dummy2);
    }
    else if(n==4) printf("\n");
    //*/
    
    //printf("before:\n");
    //print_cMatrix(&mc->accumulated_g_matsubara.matrices[n]);
    cMatrixMatrixAdditionInPlace(&mc->accumulated_g_matsubara.matrices[n], &mc->dummy2, 1.0, 1.0); // g += dummy2
    //printf("after:\n");
    //print_cMatrix(&mc->accumulated_g_matsubara.matrices[n]);
  }
  
  mc->accumulated_sign +=mc->sign;
  mc->accumulated_expOrder +=N;
  return 1;
}



int outputMeasure(MonteCarlo * mc, unsigned int nSamples) {
  int n;
  for(n=0; n<N_PTS_MAT; n++) scale_cMatrix(&mc->accumulated_g_matsubara.matrices[n],1.0/nSamples);
  
  FILE *fileOut1 = fopenSafe("green0.dat", "w", 1);
  writeToFile_cMatrixFunction(fileOut1, &mc->g0_matsubara, &mc->model);
  FILE *fileOut2 = fopenSafe("greenI.dat", "w", 1);
  
  writeToFile_cMatrixFunction(fileOut2, &mc->accumulated_g_matsubara, &mc->model);
  //writeToFile_cMatrixFunction(FILE *fileOut, cMatrixFunction * cMatFun, Model * model) {


  dMatrixFunction g_tau;
  init_dMatrixFunction(&g_tau, &mc->model);
  calculate_G0_tau(&mc->accumulated_g_matsubara,&g_tau); // warning, this function alter mc->accumulated_g_matsubara

  FILE *fileOut3 = fopenSafe("greenI_tau.dat", "w", 1);
  writeToFile_dMatrixFunction(fileOut3, &g_tau, &mc->model);


  double occupation =  -(mc->KDirac/(mc->model.beta*nSamples*mc->model.sites.n) + (mc->model.muAux - mc->model.mu) )/mc->model.U;  
  printf("nSamples = %d, sign= %f, expOrder= %f, n=%f, KDirac=%f, occupation=%f\n",
         nSamples,mc->accumulated_sign/nSamples, mc->accumulated_expOrder/nSamples, ELEM_VAL(mc->accumulated_g_tau,0,0)/nSamples, occupation, mc->density/nSamples);
  
  free_dMatrixFunction(&g_tau);
  

  
  return 1;
}


