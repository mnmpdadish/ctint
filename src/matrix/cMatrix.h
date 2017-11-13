#pragma once

#include <complex.h>
#include "../util/utilities.h"

// reminder:
// int * ptr;
// int const * ptrToConst;
// int * const constPtr;
// int const * const constPtrToConst;
// At first, I added const AFTER each *, but in the end, it is just more confusing than nothing.

// matrix matrix multiplication:
unsigned int zgemm_(char const*, char const*, unsigned int const*, unsigned int const*, unsigned int const*,
                  double complex const*, double complex const*, unsigned int const*, double complex const*, 
                  unsigned int const*, double complex const*, double complex *, unsigned int const*);
//unsigned int zswap_(unsigned int const*, double complex*, unsigned int const*, double complex*, unsigned int const*);
unsigned int zgetrf_(unsigned int const*, unsigned int const*, double complex const*, unsigned int const*,unsigned int*,unsigned int*);
unsigned int zgetri_(unsigned int const*, double complex*, unsigned int const*, unsigned int const*, double complex*, unsigned int const*, unsigned int*);

// matrix vector multiplication
unsigned int zgemv_(char const*, unsigned int const*, unsigned int const*, double complex const*, double complex const*, unsigned int const*, double complex const*, unsigned int const*, double complex const*, double complex*, unsigned int const*);

// dot product of vectors:
double complex zdotc_(unsigned int const*, double complex const*, unsigned int const*, double complex const*, unsigned int const*);

// copy arrays
unsigned int zcopy_(unsigned int const*, double complex const*, unsigned int const* , double complex*, unsigned int const*);

// add complex vectors (or matrices):
unsigned int zaxpy_(unsigned int const*, double complex const*, double complex const*, unsigned int const*, double complex*, unsigned int const*);

// multiplication by scalar
unsigned int zscal_(unsigned int const*, double complex const*, double complex*, unsigned int const*);

/*
void zscal_(const int*, const std::complex<double>*, std::complex<double>*, const int*);
void zgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*, const int*);
void zgemv_(const char*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, std::complex<double>*, const int*);
void zcopy_(const int*, const std::complex<double>*, const int* , std::complex<double>*, const int*);
void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
void zgesv_(const int*, const int*, std::complex<double>*, const int*, int*, std::complex<double>*, const int*, int*);
void zgetrf_(const int*, const int*, std::complex<double>*, const int*, int*, int*);
*/


//#define INIT_CAPACITY 4

typedef struct {
  unsigned int N;
  unsigned int capacity;
  double complex *data;
} cMatrix;

typedef cMatrix cVector;

unsigned int resize_cVector(cVector * x, unsigned int N) {
  while(N > x->capacity) x->data = realloc(x->data, (x->capacity *= 2) * sizeof(double complex));
  x->N=N;
  return 0;
}

unsigned int resize_cMatrix(cMatrix * A, unsigned int N) {
  while(N*N > A->capacity) A->data = realloc(A->data, (A->capacity *= 2) * sizeof(double complex));
  A->N=N;
  return 0;
}


unsigned int reset_cMatrix(cMatrix * A) {
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=0;
  return 0;
}


unsigned int init_cVector(cVector * x, unsigned int N) {
  x->capacity = INIT_CAPACITY;
  x->data = malloc(x->capacity * sizeof(double complex));
  x->N=N;
  resize_cVector(x,N);
  return 0;
}

unsigned int init_cMatrix(cMatrix * A, unsigned int N) {
  A->capacity = INIT_CAPACITY;
  A->data = malloc(A->capacity * sizeof(double complex));
  A->N=N;
  resize_cMatrix(A,N);
  return 0;
}

unsigned int free_cVector(cVector * x) {
  free(x->data);
  return 0;
}

unsigned int free_cMatrix(cMatrix * A) {
  free(A->data);
  return 0;
}



//Y==X
unsigned int areEqual_cVector(cVector const * X, cVector const * Y) {
  if(X->N != Y->N) return 0;
  unsigned int i;
  for(i=0; i<X->N; i++) if(!complexEqual(X->data[i], Y->data[i])) return 0;
  return 1;
}

//B==A
unsigned int areEqual_cMatrix(cMatrix const * A, cMatrix * B) {
  if(A->N != B->N) return 0;
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) if(!complexEqual(A->data[i], B->data[i])) return 0;
  return 1;
}


//Y=X
unsigned int copy_cVector(cVector const * X, cVector * Y) {
  //Y->N=X->N;
  resize_cVector(Y,X->N);
  memcpy(Y->data, X->data, X->N * sizeof(double complex));
  return 0;
}

//B=A
unsigned int copy_cMatrix(cMatrix const * A, cMatrix * B) {
  //B->N=A->N;
  resize_cMatrix(B,A->N);
  memcpy(B->data, A->data, A->N*A->N * sizeof(double complex));
  return 0;
}


#define ELEM(mtx, i, j) (mtx->data[j * mtx->N + i])
#define ELEM_VAL(mtx, i, j) (mtx.data[j * mtx.N + i])
//B=A
//where dim(B)=NxN and dim(A)=(A->N)x(A->N)
unsigned int copySub_cMatrix(cMatrix const * A, cMatrix * B, unsigned int N) {
  assert(N>=0);
  unsigned int copyLength = (A->N < N) ? A->N : N; //choose the smallest dimension between A->N and N. 
  resize_cMatrix(B,N);
  
  unsigned int i=0;	
  for(i = 0; i < copyLength; i++) {
    memcpy(&ELEM(B,0,i), &ELEM(A,0,i), sizeof(double complex)*copyLength);
  }
  return 0;
}


//C=factor*A+B
unsigned int cMatrixMatrixAddition(cMatrix const * A, cMatrix const * B, cMatrix * C, double complex factor) {
  assert(A->N == B->N);
  //resize_cMatrix(C,A->N);
  if(B!=C) copy_cMatrix(B,C);
  unsigned int N2=A->N*A->N;
  unsigned int one=1;
  zaxpy_(&N2, &factor, A->data, &one, C->data, &one);
  return 0;
}

//A=factorA*A+factorB*B
unsigned int cMatrixMatrixAdditionInPlace(cMatrix * A, cMatrix const * B, double complex factorA, double complex factorB) {
  assert(A->N == B->N);
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=factorA*A->data[i]+factorB*B->data[i];
  return 0;
}



//C=A*B
unsigned int cMatrixMatrixMultiplication(cMatrix const * A, cMatrix const * B, cMatrix * C) {
  assert(A->N == B->N);
  resize_cMatrix(C,A->N);
  unsigned int N=A->N;
  double complex one=1.0;
  double complex zero=0.0;
  char no = 'n';
  zgemm_(&no,&no,&N,&N,&N, &one, A->data, &N, B->data, &N, &zero, C->data, &N); 
  return 0;
}

/*
unsigned int cMatrixSwapCols(cMatrix * A, unsigned int row1, unsigned int row2) {
  unsigned int N=A->N;
  assert(row1 < N);
  assert(row2 < N);
  unsigned int one=1;
  dswap_(&N, &A->data[row1*A->N], &one, &A->data[row2*A->N], &one);
  return 0;
}

unsigned int cMatrixSwapRows(cMatrix * A, unsigned int col1, unsigned int col2) {
  unsigned int N=A->N;
  assert(col1 < N);
  assert(col2 < N);
  dswap_(&N, &A->data[col1], &N, &A->data[col2], &N);
  return 0;
}
*/

void swap_conjugate_complex(double complex * a, double complex * b)
{
    double complex temp = conj(*a);
    *a = conj(*b);
    *b = temp;
}

void swap_conjugate_complex_values(double complex a, double complex b)
{
    double complex temp = conj(a);
    a = conj(b);
    b = temp;
}


void conjugate_complex_values(double complex * a)
{
    *a = conj(*a);
}


void swap_complex(double complex * a, double complex * b)
{
    double complex temp = *a;
    *a = *b;
    *b = temp;
}

//A^T
unsigned int transpose_cMatrix(cMatrix * A) {
  unsigned int i,j;
  for(i=0;i<A->N;i++)
    for(j=0;j<i;j++){
      swap_complex(&ELEM(A,i,j),&ELEM(A,j,i));
    }
  return 0;
}

//A^dag
unsigned int dag_cMatrix(cMatrix * A) {
  transpose_cMatrix(A);
  unsigned int i,j;
  for(i=0;i<A->N;i++)
    for(j=0;j<A->N;j++){
      conjugate_complex_values(&ELEM(A,i,j));
    }
  return 0;
}



unsigned int print_cMatrix_real(cMatrix const * A) {
  if (!A){ 
    printf("oups.\n");
    return -1;
  }
  
  unsigned int i,j;
  double complex val;//, tol=0.00001;
  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++) {
      val = ELEM(A, i, j);
      printf("% 4.3f ", creal(val));
    }
    printf("\n");
  }
  fflush(stdout);
  return 0;
}


unsigned int print_cMatrix(cMatrix const * A) {
  if (!A){ 
    printf("oups.\n");
    return -1;
  }
  
  unsigned int i,j;
  double complex val;//, tol=0.00001;
  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++) {
      val = ELEM(A, i, j);
      printf("% 6.5f%+6.5fi ", creal(val), cimag(val));
    }
    printf("\n");
  }
  fflush(stdout);
  return 0;
}


unsigned int print_cVector(cVector const * X) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  unsigned int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.5f%+6.5fi ", creal(X->data[i]), cimag(X->data[i]));
  }
  printf("\n");

  fflush(stdout);
  return 0;
}




// A=>A^-1
void invert_cMatrix(cMatrix * A) {
  unsigned int INFO1=0;
  unsigned int INFO2=0;
  unsigned int IPIV[A->N];
  unsigned int nEntry=A->N*A->N;
  double complex WORK[nEntry]; 
  zgetrf_(&A->N, &A->N, &A->data[0], &A->N, &IPIV[0], &INFO1);
  zgetri_(&A->N, A->data, &A->N, &IPIV[0], &WORK[0], &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the cMatrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
}

// Y=A*X
void cMatrixVectorProduct(cMatrix const*A, cVector const*X, double complex const factor, cVector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resize_cVector(Y,N);
  //double complex one=1.0;
  double complex zero=0.0;
  char no = 'n';
  unsigned int inc=1;
  zgemv_(&no, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}

void cVectorMatrixProduct(cVector const*X, cMatrix const*A, double complex const factor, cVector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resize_cVector(Y,N);
  //double complex one=1.0;
  double complex zero=0.0;
  char yes = 't';
  unsigned int inc=1;
  zgemv_(&yes, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}


// return X.Y
double complex zScalarProduct(cVector const*X, cVector const*Y) {
  assert(X->N==Y->N);
  unsigned int inc=1;
  return zdotc_(&X->N, X->data, &inc, Y->data, &inc);
}

void scale_cVector(cVector *X, double complex const scal){
  unsigned int inc=1;
  zscal_(&X->N, &scal, X->data, &inc);
}

void scale_cMatrix(cMatrix *A, double complex const scal){
  unsigned int inc=1;
  unsigned int N=A->N*A->N;
  zscal_(&N, &scal, A->data, &inc);
}





