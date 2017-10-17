#pragma once

#include <stdio.h>
#include <stdlib.h>

#include "array.h"
//#include "../test.h"



int test_Array_double(int verbose)
{
  int i, Nerror=0;

  if(verbose) printf("Testing double array:\n");

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
  
  Array_double arr2;
  Array_double_init(&arr2);
  
  Array_double_copy(&arr, &arr2);
  //double pop2 = Array_double_pop(&arr2);
  Nerror+= !Array_double_areEqual(&arr, &arr2);

  if(verbose) {
    printf("copied:\n");
    for (i = 0; i < arr2.size; i++) printf("%4.2f ", arr2.data[i]);
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
  
  Array_int arr2;
  Array_int_init(&arr2);
  
  Array_int_copy(&arr, &arr2);
  //double pop2 = Array_int_pop(&arr2);
  Nerror+= !Array_int_areEqual(&arr, &arr2);

  if(verbose) {
    printf("copied:\n");
    for (i = 0; i < arr2.size; i++) printf("%d ", arr2.data[i]);
    printf("\n");
  }
  
  Array_int_free(&arr);
  
  return Nerror;
}



/*

int test_Array_void(int verbose)
{
  int i, Nerror=0;

  if(verbose) printf("\nTesting void array:\n");

  Array_void arr;
  Array_void_init(&arr, sizeof(double complex));
  Array_void_resize(&arr,4, sizeof(double complex));
  //double complex *p = &arr.data[0];
  for (i = 0; i < arr.size; i++) {
    arr.data[i]= i*1.5 + 64.0*I /i;
    printf("%f %f  ", real(arr.data[i]), imag(arr.data[i]) );
  }
  
    
  return Nerror;
}

*/

