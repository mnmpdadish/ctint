#include <stdio.h>
#include <stdlib.h>

#include "array.h"

int main(void)
{
  int i;


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
  Array_double_push(&arr, 4);
  Array_double_push(&arr, 13);
  
  for (i = 0; i < arr.size; i++)
    printf("%4.2f ", arr.data[i]);
  printf("\n");

  printf("poped! %4.2f\n",Array_double_pop(&arr));
   
  for (i = 0; i < arr.size ; i++)
    printf("%4.2f ", arr.data[i]);
  printf("\n");
  
  Array_double_free(&arr);


  printf("\nTesting integer array:\n");


  Array_int arrInt;
  Array_int_init(&arrInt);

  //Array_double_reserve(&v,32);

  Array_int_push(&arrInt, 34.8);
  Array_int_push(&arrInt, 2);
  Array_int_push(&arrInt, 4);
  Array_int_push(&arrInt, -34);
  
  for (i = 0; i < arrInt.size; i++)
    printf("%d ", arrInt.data[i]);
  printf("\n");

  printf("poped! %d\n",Array_int_pop(&arrInt));
   
  for (i = 0; i < arrInt.size ; i++)
    printf("%d ", arrInt.data[i]);
  printf("\n");
  

  Array_int_free(&arrInt);

}
