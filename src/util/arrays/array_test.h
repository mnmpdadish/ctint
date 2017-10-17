#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "array.h"
//#include "../test.h"



int test_Array_double(int verbose)
{
  int i, Nerror=0;

  printf("\nTesting double array:\n");

  Array_double arr;
  Array_double_init(&arr);

  //Array_double_reserve(&v,32);

  Array_double_push(&arr, 0.45);
  Array_double_push(&arr, 34);
  Array_double_push(&arr, 1e3);
  Array_double_push(&arr, 3.45);
  Array_double_push(&arr, 534);
  Array_double_push(&arr, 16e3);
  Array_double_push(&arr, 5);
  Array_double_push(&arr, 4.45);
  Array_double_push(&arr, 34);
  Array_double_push(&arr, 6);
  Array_double_push(&arr, 0.65);
  
  if(verbose) { 
    for (i = 0; i < arr.size; i++) printf("%4.2f ", arr.data[i]);
    printf("\n");
  }

  double pop = Array_double_pop(&arr);
  
   
  double solution[]={0.45, 34, 1e3, 3.45, 534, 16e3, 5, 4.45, 34, 6};
  
  for (i = 0; i < arr.size ; i++) if(!doubleEqual(arr.data[i],solution[i])) Nerror++;
  
  if(verbose) {
    printf("poped! %4.2f\n",pop);
    for (i = 0; i < arr.size; i++) printf("%4.2f ", arr.data[i]);
    printf("\n");
  }
  
  Array_double_free(&arr);
  return Nerror;
}


int test_Array_int(int verbose)
{
  int i, Nerror=0;


  if(verbose) printf("\nTesting integer array:\n");


  Array_int arr;
  Array_int_init(&arr);

  //Array_double_reserve(&v,32);

  Array_int_push(&arr, 38);
  Array_int_push(&arr, 2);
  Array_int_push(&arr, 4);
  Array_int_push(&arr, -34);
  

  if(verbose) { 
    for (i = 0; i < arr.size; i++) printf("%d ", arr.data[i]);
    printf("\n");
  }

  int pop = Array_int_pop(&arr);
  
  int solution[]={38,2,4,-34}; 
  for (i = 0; i < arr.size ; i++) if(arr.data[i]!=solution[i]) Nerror++;
  
  if(verbose) {
    printf("poped! %d\n",pop);
    for (i = 0; i < arr.size; i++) printf("%d ", arr.data[i]);
    printf("\n");
  }
  
  
  Array_int_free(&arr);
  
  return Nerror;
}


/*
int main() {
  int verbose=1;
  int Nfail=0;
  Nfail+= passOrFail("test_Array_double",   test_Array_double(verbose));
  Nfail+= passOrFail("test_Array_int",      test_Array_int(verbose));
  
  if(Nfail==0) printf("%s100%% of the test PASSED%s\n\n",COLORGREEN,COLORNORMAL);
  else printf("%sOh no. ABORT!%s\n\n",COLORRED,COLORNORMAL);
  
  return 0;
}
*/
