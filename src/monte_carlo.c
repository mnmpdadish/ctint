
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
  
  int i, N=50000;
  for(i=0;i<N;i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
  }
  //printf("\n\n\nnormal:\n");
  //Print_MonteCarlo(&mc);
  //CleanUpdate(&mc);
  //printf("clean:\n");
  //Print_MonteCarlo(&mc);
  
    
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


