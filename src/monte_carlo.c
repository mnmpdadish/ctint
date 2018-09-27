
#include <unistd.h>
#include <time.h>
#include "monte_carlo.h"

void do_update(MonteCarlo * mc) {
  double probFlip = 0.3;
  if(urng()<probFlip) {
    //printf("\n--------------->spin-flip\n");
    FlipVertex(mc);
  }
  else {
    //Print_MonteCarlo(mc);
    if(urng()<0.5) InsertVertex(mc); 
    else RemoveVertex(mc);
    //if(mc->vertices.N !=0) CleanUpdate(mc);
    
  }
  //Print_MonteCarlo(mc);
}

double timeSec(struct timespec start, struct timespec end){
  return (double) (end.tv_sec - start.tv_sec) + (double) (end.tv_nsec - start.tv_nsec) / 1000000000.0;
}

int main(int argc, char *argv[]) {

  struct timespec start, end;
  //usleep(2304040);
  

  if(argc!=3) {
    printf("program usage example:\n$ mc dmft.model 9\n\n");
    printf("This specific example requires the dmft.model file must be found in the current path.\n");
    printf("The 9 specify the iteration number. The hyb9.dat and params9 files must be found in the current path.\n\n");
    printf("If the iteration number is 0, there is no need for an hyb0.dat file.\n");
    printf("Only the self-consistency will be processed (with self=0) and hyb1.dat will be written.\n");
    exit(1);
  }
  
  unsigned long int iteration = atol(argv[2]);
  char * modelFileName = argv[1];
  
  char hybFileName[256];
  char paramsFileName[256];
  
  sprintf(hybFileName, "hyb%lu.dat", iteration); // puts string into buffer
  sprintf(paramsFileName, "params%lu", iteration); // puts string into buffer
  
  FILE * fileHyb;
  if(iteration>0) fileHyb = fopenSafe(hybFileName,  "rt",1);
  else if (iteration==0) fileHyb = NULL;
  else {
    printf("iteration number must be positive or 0\n");
    exit(1);
  }
  
  
  FILE * fileModel  = fopenSafe(modelFileName, "rt",1);
  FILE * fileParams = fopenSafe(paramsFileName, "rt",1);
  
  Model model;
  read_Model(fileModel, fileParams, &model, 0);
  MonteCarlo mc;
  init_MonteCarlo(fileHyb, &mc, &model);
  
  srand(mc.model.seed);
  int update_i;
  int nSamples=0;
  clock_gettime(CLOCK_MONOTONIC, &start);
  
  if(iteration>0){
    
    // Thermalization:    
    for(update_i=1; update_i<mc.model.nThermUpdates;  update_i++) {
    //for(update_i=1; update_i<20;  update_i++) {
      do_update(&mc);
      if(update_i % mc.model.cleanUpdate_i ==0) if(mc.vertices.N !=0) CleanUpdate(&mc);
      //printf("%d %d\n",update_i,mc.vertices.N); fflush(stdout);
      //if(update_i % 100 ==0) printf("%d %d\n",update_i,mc.vertices.N); fflush(stdout);
    }
    //exit(1);
    //printf("thermalization time=%f\n", (double)(clock() - start) / CLOCKS_PER_SEC); fflush(stdout);
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("thermalization time=%f\n", timeSec(start,end)); fflush(stdout);

    // Measurements:    
    clock_gettime(CLOCK_MONOTONIC, &start);
    for(update_i=1; update_i < mc.model.nUpdates+1;  update_i++) {
        
      do_update(&mc);
      if(update_i % mc.model.cleanUpdate_i ==0) {
        //printf("clean update ");
        //Print_MonteCarlo(&mc);
        if(mc.vertices.N !=0) CleanUpdate(&mc);
        //Print_MonteCarlo(&mc);
        //exit(1);
        //printf("done.\n");
        printf("%d",update_i); fflush(stdout);
        printf(".  sign=% 2.0f   order=%d   ", mc.sign, mc.vertices.N); fflush(stdout);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("  time=%f\n", timeSec(start,end)); fflush(stdout);
      }
      
      if(update_i % mc.model.measure_i==0) {
        printf("%d",update_i); fflush(stdout);
        measure(&mc);
        nSamples++;
        printf(".  sign=% 2.0f   order=%d   \n", mc.sign, mc.vertices.N); fflush(stdout);
        clock_gettime(CLOCK_MONOTONIC, &end);
        printf("  time=%f\n", timeSec(start,end)); fflush(stdout);
      }
      //printf("%d %d\n",update_i,mc.vertices.N); fflush(stdout);
    }
  }
  
  outputMeasure(&mc, nSamples, iteration);
  clock_gettime(CLOCK_MONOTONIC, &end);
  printf("total measurement time=%f\n", timeSec(start,end)); fflush(stdout);
  free_MonteCarlo(&mc);//vertices.N=0;
  return 0;
}


