
#include <stdio.h>
#include <stdlib.h>

#define SIZE 290000000

double nGlobal[SIZE];  // commented because with large size = Killed (wtf does that mean)

double createArrayRetunElement(int num) {
  double n[num];
  return n[num-2];
}

/*double createArrayRetunElement2(int num) {
  return nGlobal[num-2];
}*/

double createArrayRetunElement3(int num, double * array) {
  return array[num-2];
}


int main() {
  int i, N=SIZE;
  double value=0.;
  double total=0.;
  //double nLocal[SIZE]; // commented because with large size = segmentation fault
  double * nPtr;
  nPtr = (double *) malloc(SIZE * sizeof (double));
  if (nPtr) printf("Allocated %zu mbytes from %p to %p\n", SIZE*sizeof (double)/1024/1024, nPtr, nPtr + SIZE*sizeof (double));
  else {printf("Failed to allocated %zu mbytes\n", SIZE*sizeof(double)/1024/1024); return 1;}
  
  for(i=0;i<N;i++){
    //value = createArrayRetunElement(i);
    //value = createArrayRetunElement2(i);
    //value = createArrayRetunElement3(i, nLocal);
    value = createArrayRetunElement3(i, nPtr);
    total += value; // just to be sure that value is used, so optimize enter the function.
  }
  printf("\n tot = %f", total);
  return 0;
}

