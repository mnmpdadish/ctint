#pragma once 

#define ARRAY_INIT_CAPACITY 32

typedef struct {
    double *data;
    unsigned int capacity;
    unsigned int size;
} Array_double;

typedef struct {
    int *data;
    unsigned int capacity;
    unsigned int size;
} Array_int;


void Array_double_init(Array_double *);
void Array_double_reserve(Array_double *, int);
void Array_double_resize(Array_double *, int);
void Array_double_push(Array_double *, double);
double Array_double_pop(Array_double *);
void Array_double_free(Array_double *);

void Array_int_init(Array_int *);
void Array_int_reserve(Array_int *, int);
void Array_int_resize(Array_int *, int);
void Array_int_push(Array_int *, int);
int Array_int_pop(Array_int *);
void Array_int_free(Array_int *);

