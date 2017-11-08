
#include <unistd.h>
#include "monte_carlo.h"

int main(int argc, char *argv[]) {

  if(argc!=3) {
    printf("program usage example:\n$ mc dmft.model 9\n\n");
    printf("This specific example requires the dmft.model file must be found in the current path.\n");
    printf("The 9 specify the iteration number. The hyb9.dat and params9 files must be found in the current path.\n");
    exit(1);
  }
  
  unsigned long int iteration = atol(argv[2]);
  char * modelFileName = argv[1];
  
  char hybFileName[256];
  char paramsFileName[256];
  
  sprintf(hybFileName, "hyb%lu.dat", iteration); // puts string into buffer
  sprintf(paramsFileName, "params%lu", iteration); // puts string into buffer
  
  FILE * fileHyb    = fopenSafe(hybFileName,  "rt",1);
  FILE * fileModel  = fopenSafe(modelFileName, "rt",1);
  FILE * fileParams = fopenSafe(paramsFileName, "rt",1);
  
  Model model;
  read_Model(fileModel, fileParams, &model);
  MonteCarlo mc;
  init_MonteCarlo(fileHyb, &mc, &model);
  
  unsigned int seed = 10000061;
  srand(seed);
  
  int update_i; //, termalization_i = 10000, measure_i=10000, cleanUpdate_i=501;
  int nSamples=0;
  for(update_i=1; update_i<mc.model.nThermUpdates;  update_i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
    if(update_i % mc.model.cleanUpdate_i ==0){
      if(mc.vertices.N !=0) CleanUpdate(&mc);
    }
  }
  
  for(update_i=1; update_i < mc.model.nUpdates+1;  update_i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
    if(update_i % mc.model.measure_i==0) {
      //printf("%d",update_i); fflush(stdout);
      measure(&mc);
      nSamples++;
      //printf(".  sign=% 2.0f   order=%d   \n", mc.sign, mc.vertices.N); fflush(stdout);
    }
    if(update_i % mc.model.cleanUpdate_i ==0){
      if(mc.vertices.N !=0) CleanUpdate(&mc);
      printf("%d",update_i); fflush(stdout);
      printf(".  sign=% 2.0f   order=%d   \n", mc.sign, mc.vertices.N); fflush(stdout);
    }
  }
  
  outputMeasure(&mc, nSamples, iteration);
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


