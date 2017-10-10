#pragma once

#include "utilities.h"
#include "stdio.h"



typedef struct {
  double tau;
  unsigned int site;
  int auxSpin;  // 1=up, 0=down.
} Vertex;

typedef struct {
  int N;
  Vertex m_vertex[dataBufferSize1];
} MarkovChain;


int PrintMarkovChain(MarkovChain * markovChain){
  int i, N=markovChain->N;
  printf("\n"); for(i=0;i<N;i++) printf("%4.2f ",markovChain->m_vertex[i].tau); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",markovChain->m_vertex[i].site); 
  printf("\n"); for(i=0;i<N;i++) printf("%d    ",markovChain->m_vertex[i].auxSpin); 
  printf("\n");  
  return 0;
}

// add vertex at the end of the Markov chain:
int InsertVertex(MarkovChain * markovChain) { 
  assert(markovChain->N < dataBufferSize1);
  int N=markovChain->N;

  unsigned int nSites = 4;
  markovChain->m_vertex[N].tau = (double)rand()/(double)(RAND_MAX);
  markovChain->m_vertex[N].site= rand()%nSites;
  markovChain->m_vertex[N].auxSpin = rand()%2;
  
  markovChain->N++;
  
  return 0;
}

int RemoveVertex(MarkovChain * markovChain, unsigned int position) { 
  assert(markovChain->N >= position);

  markovChain->N--;  // the vertex is not deleted, it is just forgotten
  int N = markovChain->N;

  markovChain->m_vertex[position].tau = markovChain->m_vertex[N].tau;
  markovChain->m_vertex[position].site= markovChain->m_vertex[N].site;
  markovChain->m_vertex[position].auxSpin = markovChain->m_vertex[N].auxSpin;
  
  return 0;
}























