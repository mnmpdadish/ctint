#include <stdio.h>
#include <stdlib.h>

#include "array.h"

void Array_double_init(Array_double *arr)
{
    arr->capacity = ARRAY_INIT_CAPACITY;
    arr->size = 0;
    arr->data = malloc(sizeof(double) * arr->capacity);
}

void Array_double_reserve(Array_double *arr, int capacity)
{
  #ifdef DEBUG_ON
  printf("Array_double_resize: %d to %d\n", arr->capacity, capacity);
  #endif
  
  if(capacity > arr->capacity) { //only increase size of the array over time, if needed.
  
    void *data = realloc(arr->data, sizeof(double) * capacity);
    if (data) {
      arr->data = data;
      arr->capacity = capacity;
    }
    else {
      printf("Attempt to allocate memory failed (out of memory? who knows).\n");
      exit(1);
    }
    
  }
}

void Array_double_resize(Array_double *arr, int size)
{
    Array_double_reserve(arr, size);
    arr->size=size;
}

void Array_double_push(Array_double *arr, double item)
{
    if (arr->capacity == arr->size)
        Array_double_reserve(arr, arr->capacity * 2);
    arr->data[arr->size++] = item;
}


double Array_double_pop(Array_double *arr)
{
    double value = arr->data[arr->size-1];
    arr->size--;
    return value;
}

void Array_double_free(Array_double *arr)
{
    free(arr->data);
}






void Array_int_init(Array_int *arr)
{
    arr->capacity = ARRAY_INIT_CAPACITY;
    arr->size = 0;
    arr->data = malloc(sizeof(double) * arr->capacity);
}

void Array_int_reserve(Array_int *arr, int capacity)
{
  #ifdef DEBUG_ON
  printf("Array_int_resize: %d to %d\n", arr->capacity, capacity);
  #endif
  
  if(capacity > arr->capacity) { //only increase size of the array over time, if needed.
  
    void *data = realloc(arr->data, sizeof(double) * capacity);
    if (data) {
      arr->data = data;
      arr->capacity = capacity;
    }
    else {
      printf("Attempt to allocate memory failed (out of memory? who knows).\n");
      exit(1);
    }
    
  }
}

void Array_int_resize(Array_int *arr, int size)
{
    Array_int_reserve(arr, size);
    arr->size=size;
}

void Array_int_push(Array_int *arr, int item)
{
    if (arr->capacity == arr->size)
        Array_int_reserve(arr, arr->capacity * 2);
    arr->data[arr->size++] = item;
}


int Array_int_pop(Array_int *arr)
{
    int value = arr->data[arr->size-1];
    arr->size--;
    return value;
}

void Array_int_free(Array_int *arr)
{
    free(arr->data);
}
