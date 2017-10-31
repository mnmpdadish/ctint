
#include <unistd.h>
#include "monte_carlo.h"

int main() {

  char fileName[]="testInputFiles/plaquette2x2.model";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading model from %s:\n", fileName);
  
  Model model;
  read_Model(file,&model);
  MonteCarlo mc;
  init_MonteCarlo(&mc, &model);
  
  unsigned int seed = 1000000;
  srand(seed);
  
  int i, N=100000;
  int cleanUpdate=10000;
  for(i=0;i<N;i++) {
    double pRemove = 1.0/25.0*mc.vertices.N * (double)rand()/(double)(RAND_MAX);
    if (pRemove > 0.5) {
      int pos=rand()%mc.vertices.N;
      //printf("\n%d",pos);
      RemoveVertex(&mc, pos);
    }
    else { 
      InsertVertex(&mc);
    }
    if(i%cleanUpdate==0) Print_MonteCarlo(&mc);
    //usleep(10000);  
  }
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


