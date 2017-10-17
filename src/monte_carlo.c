
#include <unistd.h>
#include "monte_carlo.h"

int main() {

  unsigned int seed = 1000000;
  srand(seed);
  VerticesArray vertices;
  initVerticesChain(&vertices);//vertices.N=0;
  
  int i, N=100000;
  int cleanUpdate=10000;
  for(i=0;i<N;i++) {
    double pRemove = 1.0/25.0*vertices.N * (double)rand()/(double)(RAND_MAX);
    if (pRemove > 0.5) {
      int pos=rand()%vertices.N;
      //printf("\n%d",pos);
      RemoveVertex(&vertices, pos);
    }
    else { 
      InsertVertex(&vertices);
    }
    if(i%cleanUpdate==0) PrintVerticesArray(&vertices);
    //usleep(10000);  
  }
  freeVerticesChain(&vertices);//vertices.N=0;
  return 0;
}


