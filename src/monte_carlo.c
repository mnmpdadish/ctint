
#include <unistd.h>
#include "monte_carlo.h"

int main() {

  char fileName[]="testInputFiles/dmftAway.model";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found\n", fileName); exit(1);}
  printf("\nreading model from %s:\n", fileName);
  
  Model model;
  read_Model(file,&model);
  MonteCarlo mc;
  init_MonteCarlo(&mc, &model);
  
  unsigned int seed = 1000001;
  srand(seed);
  
  int update_i, N_updates=10000000, termalization_i = 10000, measure_i=10000, cleanUpdate_i=501;
  int nSamples=0;
  for(update_i=1; update_i<termalization_i;  update_i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
    if(update_i % cleanUpdate_i ==0){
      if(mc.vertices.N !=0) CleanUpdate(&mc);
    }
  }
  
  for(update_i=1; update_i<N_updates+1;  update_i++) {
    if(urng()<0.5) InsertVertex(&mc);
    else RemoveVertex(&mc);
    if(update_i % measure_i==0) {
      printf("%d",update_i); fflush(stdout);
      measure(&mc);
      nSamples++;
      printf(".  sign=% 2.0f   order=%d   \n", mc.sign, mc.vertices.N); fflush(stdout);
    }
    if(update_i % cleanUpdate_i ==0){
      if(mc.vertices.N !=0) CleanUpdate(&mc);
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
  
  int n;
  for(n=0; n<N_PTS_MAT; n++) scale_cMatrix(&mc.accumulated_g_matsubara.matrices[n],1.0/nSamples);
  
  FILE *fileOut1 = fopen("green0.dat","w");
  writeToFile_cMatrixFunction(fileOut1, &mc.g0_matsubara, &mc.model);
  FILE *fileOut2 = fopen("greenI.dat","w");
  writeToFile_cMatrixFunction(fileOut2, &mc.accumulated_g_matsubara, &mc.model);
  printf("nSamples = %d, sign= %f, expOrder= %f\n",nSamples,mc.accumulated_sign/nSamples, mc.accumulated_expOrder/nSamples);
  //writeToFile_cMatrixFunction(FILE *fileOut, cMatrixFunction * cMatFun, Model * model) {
  
    
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


