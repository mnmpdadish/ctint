#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "precalculateG0/calculateG0.h"

#define VERTICES_BASIC_CAPACITY 256

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
  dMatrixFunction g0_tau;
  //dMatrixFunction g0_tau_down;
  dMatrix M_up;   // these two matrices are THE matrices
  dMatrix M_down; // same notation as Gull. 
  dMatrix Mdummy_up;   
  dMatrix Mdummy_down; 
  dVector Q, R;
  dVector Qtilde_up, Rtilde_up; //vectors to add
  dVector Qtilde_down, Rtilde_down; 
  double S_up,   Stilde_up; //corner to add
  double S_down, Stilde_down; 
  unsigned int nSites;
  Model model;
  double sign;
} MonteCarlo;



void init_MonteCarlo(MonteCarlo * mc, Model * model) {
  init_dMatrixFunction(&mc->g0_tau, model);
  
  cMatrixFunction g0_matsubara;
  init_cMatrixFunction(&g0_matsubara, model);
  calculate_G0_matsubara(&g0_matsubara, model);
  calculate_G0_tau(&g0_matsubara,&mc->g0_tau);
  free_cMatrixFunction(&g0_matsubara);
  
  init_dMatrix(&mc->M_up,0);
  init_dMatrix(&mc->M_down,0);
  init_dMatrix(&mc->Mdummy_up,0);
  init_dMatrix(&mc->Mdummy_down,0);
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
}

void free_MonteCarlo(MonteCarlo * mc) {
  free(mc->vertices.m_vertex);
  free_dMatrixFunction(&mc->g0_tau);
  //free_dMatrixFunction(&mc->g0_tau_down);
  free_dMatrix(&mc->M_up);
  free_dMatrix(&mc->M_down);
  free_dMatrix(&mc->Mdummy_up);
  free_dMatrix(&mc->Mdummy_down);
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
  print_dMatrix(&mc->M_up);
  printf("\nM_down:\n");
  print_dMatrix(&mc->M_down);
  printf("\n");
}

double auxUp(Vertex const * vertex, MonteCarlo *mc)   { return vertex->auxSpin ? mc->model.delta : -1. - mc->model.delta;}
double auxDown(Vertex const * vertex, MonteCarlo *mc) { return vertex->auxSpin ? -1. - mc->model.delta : mc->model.delta;}
double urng() {return (double)rand()/(double)(RAND_MAX);}
double irng(unsigned int N) {return rand()%N;}


double green0(Vertex const * vertexI, Vertex const * vertexJ, MonteCarlo *mc) { 
  double diff_tau = vertexI->tau - vertexJ->tau; 

  double aps = 1.;
  if (diff_tau > .0) {
    diff_tau -= mc->model.beta;
    aps = -1.;
  }
  double ntau = fabs(diff_tau) * (N_PTS -1) /mc->model.beta;
  unsigned int ntau0= (unsigned int) ntau;
  
  double val_n  = ELEM_VAL(mc->g0_tau.matrices[ntau0], vertexI->site, vertexJ->site);
  double val_n1 = ELEM_VAL(mc->g0_tau.matrices[ntau0+1], vertexI->site, vertexJ->site);
  
  return aps*((1. - (ntau - ntau0))*val_n + (ntau - ntau0)*val_n1);
}

	

// add vertex at the end of the verticies:
int InsertVertex(MonteCarlo * mc) { 
  
  Vertex newVertex;
  newVertex.tau = mc->model.beta*urng();
  newVertex.site= irng(mc->nSites);
  newVertex.auxSpin = irng(2);
  
  mc->S_up = green0(&newVertex, &newVertex, mc) + auxUp(&newVertex, mc);
  mc->S_down = green0(&newVertex, &newVertex, mc) + auxDown(&newVertex, mc);
  printf("S_up=%f\n",mc->S_up);
  printf("S_down=%f\n",mc->S_down);
  
  if(mc->vertices.N>0){
    assert(mc->M_up.N == mc->vertices.N);
    assert(mc->M_down.N == mc->vertices.N);
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
			
    dMatrixVectorProduct(&mc->M_up, &mc->Q, 1.0, &mc->Qtilde_up);
    dMatrixVectorProduct(&mc->M_down, &mc->Q, 1.0, &mc->Qtilde_down);
    
    mc->Stilde_up = 1.0/( mc->S_up - dScalarProduct(&mc->R, &mc->Qtilde_up) );
    mc->Stilde_down = 1.0/( mc->S_down - dScalarProduct(&mc->R, &mc->Qtilde_down) );
    
    double pAcc = -2.*mc->model.sites.n*mc->model.beta*mc->model.auxU/ ( ((double)mc->vertices.N)*mc->Stilde_up*mc->Stilde_down );
    printf("pAcc=%f\n",pAcc);
    if(urng() < fabs(pAcc)) {
			
      dVectorMatrixProduct(&mc->R, &mc->M_up, -mc->Stilde_up, &mc->Rtilde_up);
      dVectorMatrixProduct(&mc->R, &mc->M_down, -mc->Stilde_down, &mc->Rtilde_down);
      
      scale_dVector(&mc->Qtilde_up, -mc->Stilde_up);
      scale_dVector(&mc->Qtilde_down, -mc->Stilde_down);
      
      dAddRowColToInverse(&mc->M_up,   &mc->Rtilde_up,   &mc->Qtilde_up,   mc->Stilde_up,   &mc->Mdummy_up);
      dAddRowColToInverse(&mc->M_down, &mc->Rtilde_down, &mc->Qtilde_down, mc->Stilde_down, &mc->Mdummy_down);
      
      copy_dMatrix(&mc->Mdummy_up,&mc->M_up);
      copy_dMatrix(&mc->Mdummy_down,&mc->M_down);

      //dMatrix swapM = &mc->M_up;
      //mc->M_up;
      
    }
    else return 0; //skip newVertex addition
  }
  else{
    double pAcc = -2.*mc->model.sites.n*mc->model.beta*mc->model.auxU*mc->S_up*mc->S_down;
    printf("%f\n",pAcc);
    if(urng() < fabs(pAcc)) {
      if(pAcc < .0) mc->sign *= -1;
      resize_dMatrix(&mc->M_up,1);
      resize_dMatrix(&mc->M_down,1);
      
      mc->M_up.data[0]   = 1./mc->S_up;
      mc->M_down.data[0] = 1./mc->S_down;
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




// remove vertex at a position (swap position with last vertex):
void RemoveVertex(MonteCarlo * mc, unsigned int position) { 
  assert(mc->vertices.N >= position);

  mc->vertices.N--;  // the vertex is not deleted, it is just forgotten
  int N = mc->vertices.N;

  mc->vertices.m_vertex[position].tau     = mc->vertices.m_vertex[N].tau;
  mc->vertices.m_vertex[position].site    = mc->vertices.m_vertex[N].site;
  mc->vertices.m_vertex[position].auxSpin = mc->vertices.m_vertex[N].auxSpin;
}







int CleanUpdate(MonteCarlo * mc) { 
  assert(mc->M_up.N == mc->vertices.N);
  assert(mc->M_down.N == mc->vertices.N);
  unsigned int i,j;
  for(i = 0; i < mc->vertices.N; i++) {
    Vertex vertexI = mc->vertices.m_vertex[i];
    for(j = 0; j < mc->vertices.N; j++) {
      Vertex vertexJ = mc->vertices.m_vertex[j];
      ELEM_VAL(mc->M_up,i,j) = green0(&vertexI, &vertexJ, mc);
      ELEM_VAL(mc->M_down,i,j) = green0(&vertexI, &vertexJ, mc);
    }
    ELEM_VAL(mc->M_up,i,i) += auxUp(&vertexI, mc);
    ELEM_VAL(mc->M_down,i,i) += auxDown(&vertexI, mc);
  }
  invert_dMatrix(&mc->M_up);
  invert_dMatrix(&mc->M_down);
  return 1;
}


