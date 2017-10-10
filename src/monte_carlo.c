
#include <unistd.h>
#include "monte_carlo.h"

int main() {

  unsigned int seed = 1000000;
  srand(seed);
  MarkovChain markovChain;
  markovChain.N=0;
  

  int i, N=100000;
  int cleanUpdate=10000;
  for(i=0;i<N;i++) {
    double pRemove = 1.0/25.0*markovChain.N * (double)rand()/(double)(RAND_MAX);
    if (pRemove > 0.5) {
      int pos=rand()%markovChain.N;
      //printf("\n%d",pos);
      RemoveVertex(&markovChain, pos);
    }
    else { 
      InsertVertex(&markovChain);
    }
    if(i%cleanUpdate==0) PrintMarkovChain(&markovChain);
    //usleep(10000);  
  }
  return 0;
}


