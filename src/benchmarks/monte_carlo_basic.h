#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
//#include "util/utilities.h"

#define VERTICES_BASIC_CAPACITY 4

typedef struct {
  double tau;
  unsigned int site;
  int auxSpin;  // 1=up, 0=down.
} Vertex;

typedef struct {
  unsigned int N;
  unsigned int capacity;
  Vertex * m_vertex;
} VerticesArray;



unsigned int initVerticesChain(VerticesArray * vertices) {
  //Array_double_init(&x->buffer);
  vertices->capacity = VERTICES_BASIC_CAPACITY;
  vertices->m_vertex = malloc(vertices->capacity * sizeof(Vertex));
  vertices->N=0;
  return 0;
}


unsigned int freeVerticesChain(VerticesArray * vertices) {
  free(vertices->m_vertex);
  return 0;
}

int PrintVerticesArray(VerticesArray * vertices){
  int i, N=vertices->N;
  printf("\n"); for(i=0;i<N;i++) printf("%4.2f ",vertices->m_vertex[i].tau); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",vertices->m_vertex[i].site); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",vertices->m_vertex[i].auxSpin); 
  printf("\n");  
  return 0;
}

// add vertex at the end of the Markov chain:
int InsertVertex(VerticesArray * vertices) { 
  //assert(vertices->N < 1000);
  if(vertices->N == vertices->capacity) {
    vertices->m_vertex = realloc(vertices->m_vertex, (vertices->capacity *= 2) * sizeof(Vertex));
    //printf("Resizing... %d %d\n",vertices->N, vertices->capacity );
  }
  
  int N=vertices->N;

  unsigned int nSites = 4;
  vertices->m_vertex[N].tau = (double)rand()/(double)(RAND_MAX);
  vertices->m_vertex[N].site= rand()%nSites;
  vertices->m_vertex[N].auxSpin = rand()%2;
  
  vertices->N++;
  
  return 0;
}

int RemoveVertex(VerticesArray * vertices, unsigned int position) { 
  assert(vertices->N >= position);

  vertices->N--;  // the vertex is not deleted, it is just forgotten
  int N = vertices->N;

  vertices->m_vertex[position].tau = vertices->m_vertex[N].tau;
  vertices->m_vertex[position].site= vertices->m_vertex[N].site;
  vertices->m_vertex[position].auxSpin = vertices->m_vertex[N].auxSpin;
  
  return 0;
}









