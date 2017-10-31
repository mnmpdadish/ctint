
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
  
  int i;
  //int cleanUpdate=10;
  for(i=0;i<1000;i++) {
    InsertVertex(&mc);
/*    if(i%cleanUpdate==0){
      printf("\n\n\nnormal:\n");
      Print_MonteCarlo(&mc);
      CleanUpdate(&mc);
      printf("clean:\n");
      Print_MonteCarlo(&mc);
    }
*/  
    //usleep(10000);  
  }
  for(i=0;i<995;i++) {
    RemoveVertex(&mc);
    
    //usleep(10000);  
  }
  printf("\n\n\nnormal:\n");
  Print_MonteCarlo(&mc);
  CleanUpdate(&mc);
  printf("clean:\n");
  Print_MonteCarlo(&mc);
  
    
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


