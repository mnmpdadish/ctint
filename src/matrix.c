#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


extern int dgemm_(char*,char*,int*,int*,int*,double*,double*,int*,double*, int*, double*,double*, int*);


typedef struct {
  int N;
  int dataBufferSize; // can have a much larger buffer, in order to enable easy resizing
  double * data;
} squareMatrix;

squareMatrix * newSquareMatrix(int dataBufferSize) {
  if (dataBufferSize <= 0) return NULL;

  // allocate a squareMatrix structure
  squareMatrix * m = (squareMatrix *) malloc(sizeof(squareMatrix));
  m->N = 0;
  m->dataBufferSize = dataBufferSize;
  m->data = (double *) malloc(m->dataBufferSize*sizeof(double));
  // set all data to 0
  int i;
  for (i = 0; i < dataBufferSize; i++)
    m->data[i] = 0.0;
  return m;
}

int deleteMatrix(squareMatrix * mtx) {
  if (!mtx) return -1;
  // free mtx's data
  assert (mtx->data);
  free(mtx->data);
  // free mtx itself
  free(mtx);
  return 0;
}

int resizeMatrix(squareMatrix * mtx, int N) {
  assert(N >= 0);
  assert(N*N <= mtx->dataBufferSize);
  mtx->N=N;
  //memcpy(cp->data, mtx->data, mtx->rows * mtx->cols * sizeof(double));
  return 0;
}

int copySquareMatrix(squareMatrix * mtxIN, squareMatrix * mtxOUT) {
  assert(mtxOUT->dataBufferSize >= mtxIN->dataBufferSize);
  mtxOUT->N=mtxIN->N;
  memcpy(mtxOUT->data, mtxIN->data, mtxIN->N * mtxIN->N * sizeof(double));
  return 0;
}

int matrixMatrixMultiplication(squareMatrix * A, squareMatrix * B, squareMatrix * C) {
  assert(C->dataBufferSize >= A->dataBufferSize);
  assert(C->dataBufferSize >= B->dataBufferSize);
  assert(A->N == B->N);
  C->N=A->N;
  int N=A->N;
  double one=1.0;
  double zero=1.0;
  char yes = 't';
  dgemm_(&yes,&yes,&N,&N,&N, &one, A->data, &N, B->data, &N, &zero, C->data, &N); //we need to transpose because lapack does not define the matrices the same way we do.
  return 0;
}


/*
        zgetrf_(&dim,&dim,data_,&dim,IPIV_,&INFO);
        zgetri_(&dim,data_,&dim,IPIV_,WORK_,&nEntry,&INFO);
        if( !(INFO == 0) ) {
            printf( "The algorithm failed to invert the matrix.\n" );
            exit( 1 );
*/

#define ELEM(mtx, i, j) \
  mtx->data[i * mtx->N + j]

int printMatrix(squareMatrix * mtx) {
  if (!mtx) return -1;
  
  int i,j;
  for (i = 0; i < mtx->N; i++) {
    for (j = 0; j < mtx->N; j++) {
      printf("% 6.2f ", ELEM(mtx, i, j));
    }
    printf("\n");
  }
  return 0;
}


int testCopy() {
  squareMatrix *A, *B;
  A = newSquareMatrix(100*100);
  B = newSquareMatrix(100*100);
  
  resizeMatrix(A,3);
  
  double * p = A->data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  printf("\nA=\n"); printMatrix(A);
  printf("\nB=\n"); printMatrix(B);

  copySquareMatrix(A,B);
  printf("\ncopying\nB=\n"); printMatrix(B);
  
  deleteMatrix(A);
  deleteMatrix(B);
  return 0;
}


int testMultiply() {
  squareMatrix *A, *B, *C;
  A = newSquareMatrix(100*100);
  B = newSquareMatrix(100*100);
  C = newSquareMatrix(100*100);
  
  resizeMatrix(A,3);
  resizeMatrix(B,3);
  
  double * p = A->data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  
  p = B->data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 2.0; *p++ = 2.0;
  
  printf("\nA=\n"); printMatrix(A);
  printf("\nB=\n"); printMatrix(B);

  matrixMatrixMultiplication(A,B,C);
  printf("\nC=A*B=\n"); printMatrix(C);

  deleteMatrix(A);
  deleteMatrix(B);
  deleteMatrix(C);
  return 0;
}



int main() {
  //printf("\n---------\ntestCopy():\n");
  //testCopy();
  printf("\n---------\ntestMultiply():\n");
  testMultiply();


  return 0;
}
