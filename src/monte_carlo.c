
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
  
  int i, N=11;
  int cleanUpdate=10;
  for(i=0;i<N;i++) {
    //double pRemove = (double)rand()/(double)(RAND_MAX);
    InsertVertex(&mc);
    printf("\n\n\nnormal:\n");
    Print_MonteCarlo(&mc);
    if(i%cleanUpdate==0){
      CleanUpdate(&mc);
      printf("clean:\n");
      Print_MonteCarlo(&mc);
    }
    
    //usleep(10000);  
  }
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


