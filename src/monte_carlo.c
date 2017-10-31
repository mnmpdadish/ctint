
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
  
  int update_i, N=10000, measure_i=500;
  for(update_i=0;update_i<N;update_i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
    if(update_i % measure_i) {
      //measure(&mc);
    }
  }
  //printf("\n\n\nnormal:\n");
  //Print_MonteCarlo(&mc);
  //CleanUpdate(&mc);
  //printf("clean:\n");
  //Print_MonteCarlo(&mc);
  
    
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


