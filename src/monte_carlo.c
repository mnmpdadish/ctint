
#include <unistd.h>
#include "monte_carlo.h"

int main() {

  //char fileName[]="testInputFiles/dmftAway.model";
  FILE * file = fopenSafe("files/dmftAway.model", "rt",1);
  //if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  //printf("\nreading model from %s:\n", fileName);
  
  Model model;
  read_Model(file,&model);
  MonteCarlo mc;
  init_MonteCarlo(&mc, &model);
  
  unsigned int seed = 100000;
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
    /*
    printf("\n\n\nnormal:\n");
    Print_MonteCarlo(&mc);
    if(mc.vertices.N !=0) CleanUpdate(&mc);
    printf("clean:\n");
    Print_MonteCarlo(&mc);
    printf("sign=% 2.0f   order=%d   \n\n", mc.sign, mc.vertices.N); fflush(stdout);
    //*/
  }
  
  outputMeasure(&mc, nSamples);
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


