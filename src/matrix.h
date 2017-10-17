#pragma once

#include "util/utilities.h"
//#include "matrix.h"

#define INIT_CAPACITY 4

typedef struct {
  unsigned int N;
  unsigned int capacity;
  double *data;
} matrix;

typedef struct {
  unsigned int N;
  unsigned int capacity;
  double *data;
  //double data[DATA_BUFFER_SIZE_1];
} vector;

unsigned int resizeVector(vector * x, unsigned int N) {
  //assert(N >= 0);
  //assert(N <= DATA_BUFFER_SIZE_1);
  while(N > x->capacity) x->data = realloc(x->data, (x->capacity *= 2) * sizeof(double));
  x->N=N;
  return 0;
}

unsigned int resizeMatrix(matrix * A, unsigned int N) {
  //assert(N >= 0);
  //assert(N*N <= DATA_BUFFER_SIZE_2);
  while(N*N > A->capacity) A->data = realloc(A->data, (A->capacity *= 2) * sizeof(double));
  A->N=N;
  return 0;
}

unsigned int resetMatrix(matrix * A) {
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=0;
  return 0;
}


unsigned int initVector(vector * x, unsigned int N) {
  //Array_double_init(&x->buffer);
  x->capacity = INIT_CAPACITY;
  x->data = malloc(x->capacity * sizeof(double));
  x->N=N;
  resizeVector(x,N);
  return 0;
}

unsigned int initMatrix(matrix * A, unsigned int N) {
  //Array_double_init(&A->buffer);
  A->capacity = INIT_CAPACITY;
  A->data = malloc(A->capacity * sizeof(double));
  A->N=N;
  resizeMatrix(A,N);
  return 0;
}

unsigned int freeVector(vector * x) {
  //Array_double_free(&x->buffer);
  free(x->data);
  return 0;
}

unsigned int freeMatrix(matrix * A) {
  //Array_double_free(&A->buffer);
  free(A->data);
  return 0;
}



//Y==X
unsigned int areEqual_V(vector const * X, vector const * Y) {
  if(X->N != Y->N) return 0;
  unsigned int i;
  for(i=0; i<X->N; i++) if(!doubleEqual(X->data[i], Y->data[i])) return 0;
  return 1;
}

//B==A
unsigned int areEqual_M(matrix const * A, matrix * B) {
  if(A->N != B->N) return 0;
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) if(!doubleEqual(A->data[i], B->data[i])) return 0;
  return 1;
}


//Y=X
unsigned int copyVector(vector const * X, vector * Y) {
  //Y->N=X->N;
  resizeVector(Y,X->N);
  memcpy(Y->data, X->data, X->N * sizeof(double));
  return 0;
}

//B=A
unsigned int copyMatrix(matrix const * A, matrix * B) {
  //B->N=A->N;
  resizeMatrix(B,A->N);
  memcpy(B->data, A->data, A->N*A->N * sizeof(double));
  return 0;
}


#define ELEM(mtx, i, j) (mtx->data[j * mtx->N + i])
//B=A
//where dim(B)=NxN and dim(A)=(A->N)x(A->N)
unsigned int copySubMatrix(matrix const * A, matrix * B, unsigned int N) {
  assert(N>=0);
  unsigned int copyLength = (A->N < N) ? A->N : N; //choose the smallest dimension between A->N and N. 
  resizeMatrix(B,N);
  
  unsigned int i=0;	
  for(i = 0; i < copyLength; i++) {
    memcpy(&ELEM(B,0,i), &ELEM(A,0,i), sizeof(double)*copyLength);
  }
  return 0;
}


//C=A*B
unsigned int matrixMatrixMultiplication(matrix const * A, matrix const * B, matrix * C) {
  assert(A->N == B->N);
  resizeMatrix(C,A->N);
  unsigned int N=A->N;
  double one=1.0;
  double zero=0.0;
  char no = 'n';
  dgemm_(&no,&no,&N,&N,&N, &one, A->data, &N, B->data, &N, &zero, C->data, &N); 
  return 0;
}

unsigned int matrixSwapCols(matrix * A, unsigned int row1, unsigned int row2) {
  unsigned int N=A->N;
  assert(row1 < N);
  assert(row2 < N);
  unsigned int one=1;
  dswap_(&N, &A->data[row1*A->N], &one, &A->data[row2*A->N], &one);
  return 0;
}

unsigned int matrixSwapRows(matrix * A, unsigned int col1, unsigned int col2) {
  unsigned int N=A->N;
  assert(col1 < N);
  assert(col2 < N);
  dswap_(&N, &A->data[col1], &N, &A->data[col2], &N);
  return 0;
}

void swap_double(double * a, double * b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

void swap_ptrOfDoubles(double * a, double * b)
{
    double * temp = a;
    a = b;
    b = temp;
}


//A^T
unsigned int transposeMatrix(matrix * A) {
  unsigned int i,j;
  for(i=0;i<A->N;i++)
    for(j=0;j<i;j++){
      swap_double(&ELEM(A,i,j),&ELEM(A,j,i));
    }
  return 0;
}


unsigned int printMatrix(matrix const * A) {
  if (!A){ 
    printf("oups.\n");
    return -1;
  }
  
  unsigned int i,j;
  double val;//, tol=0.00001;
  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++) {
      val = ELEM(A, i, j);
      //if(fabs(val)>tol) 
      printf("% 6.3f ", val);
      //else printf("   .   ");
    }
    printf("\n");
  }
  fflush(stdout);
  return 0;
}


unsigned int printVector(vector const * X) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  unsigned int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


unsigned int printVectorFactor(vector const * X, double f) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  unsigned int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", f*X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


// A=>A^-1
void invertMatrix(matrix * A) {
  unsigned int INFO1=0;
  unsigned int INFO2=0;
  unsigned int IPIV[A->N];
  unsigned int nEntry=A->N*A->N;
  double WORK[nEntry]; 
  dgetrf_(&A->N, &A->N, &A->data[0], &A->N, &IPIV[0], &INFO1);
  dgetri_(&A->N, A->data, &A->N, &IPIV[0], &WORK[0], &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the matrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
}

// Y=A*X
void matrixVectorProduct(matrix const*A, vector const*X, double const factor, vector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resizeVector(Y,N);
  //double one=1.0;
  double zero=0.0;
  char no = 'n';
  unsigned int inc=1;
  dgemv_(&no, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}

void vectorMatrixProduct(vector const*X, matrix const*A, double const factor, vector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resizeVector(Y,N);
  //double one=1.0;
  double zero=0.0;
  char yes = 't';
  unsigned int inc=1;
  dgemv_(&yes, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}


// return X.Y
double scalarProduct(vector const*X, vector const*Y) {
  assert(X->N==Y->N);
  unsigned int inc=1;
  return ddot_(&X->N, X->data, &inc, Y->data, &inc);
}

void scaleVector(vector *X, double const scal){
  unsigned int inc=1;
  dscal_(&X->N, &scal, X->data, &inc);
}



// Suppose a matrix A composed of the block matrices Aij:
// A = [ A11 A12 ]
//     [ A21 A22 ]
//
// which is the inverse of D, such that:
// [ D11 D12 ]  [ A11 A12 ] = [ 1 0 ]
// [ D21 D22 ]  [ A11 A12 ]   [ 0 1 ]
//
// then we can prove directly that the Schur complement is:
// D11^(-1) = A11 - A12 A22^(-1) A21
//
// At the output of the function S = D11^-1 (the Schur complement).
// Here, dimensions are:
// dim(A) = NxN, 
// dim(A11)=(N-1)x(N-1), dim(A12)=1x(N-1), dim(A21)=(N-1)xN, dim(A22)=1x1
// Same for Dij and dim(S)=(N-1)x(N-1).
unsigned int schurComplement(matrix const*A, matrix *S) {
  unsigned int inc = 1;
  unsigned int N = A->N;
  unsigned int Nm1 = N-1;
  //S->N=Nm1;
  //assert(S->N==N-1);
  assert(Nm1>0);
  copySubMatrix(A,S,Nm1);

  double factor = -1./ELEM(A,Nm1,Nm1);
  // this next line does S = S - A12 A22^(-1) A21;
  dger_(&S->N, &S->N, &factor, &ELEM(A,0,Nm1), &inc, &ELEM(A,Nm1,0), &A->N, S->data, &S->N);
  return 0;
}



unsigned int addRowColToInverse(matrix const*A, vector const*Rtilde, vector const*Qtilde, double const sTilde, matrix *Ap1) {
  unsigned int inc = 1;
  unsigned int N = A->N;
  assert(N==Rtilde->N);
  assert(N==Qtilde->N);
  unsigned int Np1 = N+1;
  copySubMatrix(A,Ap1,Np1);

  //double one  = 1.0;
  //double minus=-1.0;
  //double zero = 0.0;
  
  double factor = 1.0/sTilde;
  // got pAccept now.
  dger_(&N, &N, &factor, Qtilde->data, &inc, Rtilde->data, &inc, Ap1->data, &Np1);
  //printf ("\nAp1=\n"); printMatrix(Ap1);
  
  //daxpy_(&N, &factor, Qtilde->data, &inc, &ELEM(Ap1,0,N), &inc);				
  //printf ("\nAp1=\n"); printMatrix(Ap1);
  
  dcopy_(&N, Rtilde->data, &inc, &ELEM(Ap1,N,0), &Np1);
  dcopy_(&N, Qtilde->data, &inc, &ELEM(Ap1,0,N), &inc);
  
  
  //printf ("\nAp1=\n"); printMatrix(Ap1);
  
  ELEM(Ap1,N,N)=sTilde;
  
  
  
  

  //double factor = -1./ELEM(A,Nm1,Nm1);
  // this next line does S = S - A12 A22^(-1) A21;
  //dger_(&S->N, &S->N, &factor, &ELEM(A,0,Nm1), &inc, &ELEM(A,Nm1,0), &A->N, S->data, &S->N);
  return 0;
}
		































