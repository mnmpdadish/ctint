#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define N 10

typedef struct {
  int * array1;
  int * array2;
} Basic;


// this is juste a basic test
int main() {
  int i;
  Basic basic;
  basic.array1 = malloc(N*sizeof(int));
  basic.array2 = malloc(2*N*sizeof(int));
  Basic basic2 = basic;
  for(i=0;i<N;i++){
    basic.array1[i] = i*i;
    basic.array2[2*i] = i;
    basic.array2[2*i+1] = i*i-i;
    printf("%d %d %d\n", basic.array1[i], basic.array2[2*i], basic.array2[2*i+1]);
  }
  
  printf("\ncopying struct...\n");
  for(i=0;i<N;i++){
    printf("%d %d %d\n", basic2.array1[i], basic2.array2[2*i], basic2.array2[2*i+1]);
  }
  
  return 0;
}

