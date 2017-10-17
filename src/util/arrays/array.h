#pragma once 

#include "../utilities.h"

#define ARRAY_INIT_CAPACITY 4


typedef struct {
    double *data;
    unsigned int capacity;
    unsigned int size;
} Array_double;


// ---------------------------------------------------------------
// -------------------------- DOUBLE -----------------------------
// ---------------------------------------------------------------


void Array_double_init(Array_double *arr)
{
    arr->capacity = ARRAY_INIT_CAPACITY;
    arr->size = 0;
    arr->data = malloc(sizeof(double) * arr->capacity);
}

void Array_double_reserve(Array_double *arr, int capacity)
{
  
  if(capacity > arr->capacity) { //only increase size of the array over time, if needed.
    #ifdef DEBUG_ON
    printf("Array_double_resize: %d to %d\n", arr->capacity, capacity);
    #endif
  
    arr->data = realloc(arr->data, sizeof(double) * capacity);
    arr->capacity = capacity;
    if (!arr->data) {
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

//Y==X
int Array_double_areEqual(Array_double const *arr1, Array_double const *arr2) {
  if(arr1->size != arr2->size) return 0;
  int i;
  for(i=0; i<arr1->size; i++) if(!doubleEqual(arr1->data[i], arr2->data[i])) return 0;
  return 1;
}

int Array_double_copy(Array_double const * A, Array_double * B) {
  Array_double_resize(B,A->size);
  memcpy(B->data, A->data, A->size * sizeof(double));
  return 0;
}



// ---------------------------------------------------------------
// -------------------------- INT --------------------------------
// ---------------------------------------------------------------


typedef struct {
    int *data;
    unsigned int capacity;
    unsigned int size;
} Array_int;


void Array_int_init(Array_int *arr)
{
    arr->capacity = ARRAY_INIT_CAPACITY;
    arr->size = 0;
    arr->data = malloc(sizeof(int) * arr->capacity);
}


void Array_int_reserve(Array_int *arr, int capacity)
{
  
  if(capacity > arr->capacity) { //only increase size of the array over time, if needed.
    #ifdef DEBUG_ON
    printf("Array_int_resize: %d to %d\n", arr->capacity, capacity);
    #endif
    
    arr->data = realloc(arr->data, sizeof(int) * capacity);
    arr->capacity = capacity;
    if (!arr->data) {
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

//Y==X
int Array_int_areEqual(Array_int const *arr1, Array_int const *arr2) {
  if(arr1->size != arr2->size) return 0;
  int i;
  for(i=0; i<arr1->size; i++) if(arr1->data[i] != arr2->data[i]) return 0;
  return 1;
}

int Array_int_copy(Array_int const * A, Array_int * B) {
  Array_int_resize(B,A->size);
  memcpy(B->data, A->data, A->size * sizeof(int));
  return 0;
}




/*
// ---------------------------------------------------------------
// -------------------------- VOID -------------------------------
// ---------------------------------------------------------------


typedef struct {
    void *data;
    unsigned int capacity;
    unsigned int size;
} Array_void;


void Array_void_init(Array_void *arr, int sizeofType)
{
    arr->capacity = ARRAY_INIT_CAPACITY;
    arr->size = 0;
    arr->data = malloc(sizeofType * arr->capacity);
}

void Array_void_reserve(Array_void *arr, int capacity, int sizeofType)
{
  
  if(capacity > arr->capacity) { //only increase size of the array over time, if needed.
    #ifdef DEBUG_ON
    printf("Array_void_resize: %d to %d\n", arr->capacity, capacity);
    #endif
  
    void *data = realloc(arr->data, sizeofType * capacity);
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

void Array_void_resize(Array_void *arr, int size, int sizeofType)
{
    Array_void_reserve(arr, size, sizeofType);
    arr->size=size;
}

void Array_void_free(Array_void *arr)
{
    free(arr->data);
}

*/
