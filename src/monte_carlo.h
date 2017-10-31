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
  dMatrixFunction g0_tau_up;
  dMatrixFunction g0_tau_down;
  dMatrix M_up;   // these two matrices are THE matrices
  dMatrix M_down; // same notation as Gull. 
  dVector Q, R, Qtilde, Rtilde; //vectors to add
  double S, Stilde; //corner to add
  unsigned int nSites;
} MonteCarlo;



void init_MonteCarlo(MonteCarlo * mc, Model * model) {
  init_dMatrixFunction(&mc->g0_tau_up, model);
  init_dMatrixFunction(&mc->g0_tau_down, model);
  
  init_dMatrix(&mc->M_up,0);
  init_dMatrix(&mc->M_down,0);
  init_dVector(&mc->Q,0);
  init_dVector(&mc->R,0);
  init_dVector(&mc->Qtilde,0);
  init_dVector(&mc->Rtilde,0);
  
  mc->nSites = model->sites.n;
  
  mc->vertices.capacity = VERTICES_BASIC_CAPACITY;
  mc->vertices.m_vertex = malloc(mc->vertices.capacity * sizeof(Vertex));
  mc->vertices.N=0;
}

void free_MonteCarlo(MonteCarlo * mc) {
  free(mc->vertices.m_vertex);
  free_dMatrixFunction(&mc->g0_tau_up);
  free_dMatrixFunction(&mc->g0_tau_down);
  free_dMatrix(&mc->M_up);
  free_dMatrix(&mc->M_down);
  free_dMatrix(&mc->Q);
  free_dMatrix(&mc->R);
  free_dMatrix(&mc->Qtilde);
  free_dMatrix(&mc->Rtilde);
}

void Print_MonteCarlo(MonteCarlo * mc){
  int i, N=mc->vertices.N;
  printf("\n"); for(i=0;i<N;i++) printf("%4.2f ",mc->vertices.m_vertex[i].tau); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].site); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",mc->vertices.m_vertex[i].auxSpin); 
  printf("\n");
  print_dMatrix(&mc->M_up);
  printf("\n");
  print_dMatrix(&mc->M_down);
  printf("\n");
}

// add vertex at the end of the verticies:
void InsertVertex(MonteCarlo * mc) { 
  if(mc->vertices.N == mc->vertices.capacity) {
    mc->vertices.m_vertex = realloc(mc->vertices.m_vertex, (mc->vertices.capacity *= 2) * sizeof(Vertex));
  }
  int N=mc->vertices.N;
  
  if(N>0){
    dMatrixVectorProduct(&mc->M_up, &mc->Q, 1.0, &mc->Qtilde);
    mc->Stilde = 1.0/( mc->S - dScalarProduct(&mc->R,&mc->Qtilde) );
    dVectorMatrixProduct(&mc->R, &mc->M_up, -mc->Stilde,&mc->Rtilde);
    scale_dVector(&mc->Qtilde,-mc->Stilde);
  }
  else{
    resize_dMatrix(&mc->M_up,N+1);
    resize_dMatrix(&mc->M_down,N+1);
    
  }
  
  mc->vertices.m_vertex[N].tau = (double)rand()/(double)(RAND_MAX);
  mc->vertices.m_vertex[N].site= rand()%mc->nSites;
  mc->vertices.m_vertex[N].auxSpin = rand()%2;
  
  mc->vertices.N++;
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

