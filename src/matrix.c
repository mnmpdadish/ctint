#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// remainder:
// int * ptr;
// int const * ptrToConst;
// int * const constPtr;
// int const * const constPtrToConst;
// At first, I added const AFTER each *, but in the end, it is just more confusing than nothing.


int dgemm_(char const*, char const*, int const*, int const*, int const*,
                  double const*, double const*, int const*, double const*, 
                  int const*, double const*, double *, int const*);
int dswap_(int const*, double*, int const*, double*, int const*);
int dgetrf_(int const*, int const*, double const*, int const*,int*,int*);
int dgetri_(int const*, double*, int const*, int const*, double*, int const*, int*);
int dgemv_(char const*, int const*, int const*, double const*, double const*, int const*, double const*, int const*, double const*, double*, int const*);
double ddot_(int const*, double const*, int const*, double const*, int const*);
int dger_(int const*, int const*, double const*, double const*, int const*, double const*, int const*, double *, int const*);
int dcopy_(int const*, double const*, int const* , double*, int const*);
int daxpy_(int const*, double const*, double const*, int const*, double*, int const*);
int dscal_(int const*, double const*, double*, int const*);


	
double fabs(double);
  
//#define DATA_BUFFER_SIZE_1 100
//#define DATA_BUFFER_SIZE_2 10000 // must be square of DATA_BUFFER_SIZE_1

typedef struct {
  int N;
  double data[DATA_BUFFER_SIZE_2];
} matrix;

typedef struct {
  int N;
  double data[DATA_BUFFER_SIZE_1];
} vector;

int resizeVector(vector * x, int N) {
  assert(N >= 0);
  assert(N <= DATA_BUFFER_SIZE_1);
  x->N=N;
  return 0;
}

int resizeMatrix(matrix * A, int N) {
  assert(N >= 0);
  assert(N*N <= DATA_BUFFER_SIZE_2);
  A->N=N;
  return 0;
}

int resetMatrix(matrix * A) {
  int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=0;
  return 0;
}


int doubleEqual(double const a, double const b) {
    return fabs(a - b) < 0.000000001;
}

//Y==X
int areEqual_V(vector const * X, vector const * Y) {
  if(X->N != Y->N) return 0;
  int i;
  for(i=0; i<X->N; i++) if(!doubleEqual(X->data[i], Y->data[i])) return 0;
  return 1;
}

//B==A
int areEqual_M(matrix const * A, matrix * B) {
  if(A->N != B->N) return 0;
  int i;
  for(i=0; i<(A->N*A->N); i++) if(!doubleEqual(A->data[i], B->data[i])) return 0;
  return 1;
}


//Y=X
int copyVector(vector const * X, vector * Y) {
  Y->N=X->N;
  memcpy(Y->data, X->data, X->N * sizeof(double));
  return 0;
}

//B=A
int copyMatrix(matrix const * A, matrix * B) {
  B->N=A->N;
  memcpy(B->data, A->data, A->N*A->N * sizeof(double));
  return 0;
}


#define ELEM(mtx, i, j) (mtx->data[j * mtx->N + i])
//B=A
//where dim(B)=NxN and dim(A)=(A->N)x(A->N)
int copySubMatrix(matrix const * A, matrix * B, int N) {
  assert(N>=0);
  int copyLength = (A->N < N) ? A->N : N; //choose the smallest dimension between A->N and N. 
  B->N=N;
  
  int i=0;	
  for(i = 0; i < copyLength; i++) {
    memcpy(&ELEM(B,0,i), &ELEM(A,0,i), sizeof(double)*copyLength);
  }
  return 0;
}


//C=A*B
int matrixMatrixMultiplication(matrix const * A, matrix const * B, matrix * C) {
  assert(A->N == B->N);
  C->N=A->N;
  int N=A->N;
  double one=1.0;
  double zero=0.0;
  char no = 'n';
  dgemm_(&no,&no,&N,&N,&N, &one, A->data, &N, B->data, &N, &zero, C->data, &N); 
  return 0;
}

int matrixSwapCols(matrix * A, int row1, int row2) {
  int N=A->N;
  assert(row1 < N);
  assert(row2 < N);
  int one=1;
  dswap_(&N, &A->data[row1*A->N], &one, &A->data[row2*A->N], &one);
  return 0;
}

int matrixSwapRows(matrix * A, int col1, int col2) {
  int N=A->N;
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
int transposeMatrix(matrix * A) {
  int i,j;
  for(i=0;i<A->N;i++)
    for(j=0;j<i;j++){
      swap_double(&ELEM(A,i,j),&ELEM(A,j,i));
    }
  return 0;
}


int printMatrix(matrix const * A) {
  if (!A){ 
    printf("oups.\n");
    return -1;
  }
  
  int i,j;
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


int printVector(vector const * X) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


int printVectorFactor(vector const * X, double f) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", f*X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


// A=>A^-1
void invertMatrix(matrix * A) {
  int INFO1=0;
  int INFO2=0;
  int IPIV[A->N];
  int nEntry=A->N*A->N;
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
  int N=A->N;
  assert(N == X->N);
  assert(N == Y->N);
  //double one=1.0;
  double zero=0.0;
  char no = 'n';
  int inc=1;
  dgemv_(&no, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}

void vectorMatrixProduct(vector const*X, matrix const*A, double const factor, vector *Y) {
  int N=A->N;
  assert(N == X->N);
  assert(N == Y->N);
  //double one=1.0;
  double zero=0.0;
  char yes = 't';
  int inc=1;
  dgemv_(&yes, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}


// return X.Y
double scalarProduct(vector const*X, vector const*Y) {
  assert(X->N==Y->N);
  int inc=1;
  return ddot_(&X->N, X->data, &inc, Y->data, &inc);
}

void scaleVector(vector *X, double const scal){
  int inc=1;
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
int schurComplement(matrix const*A, matrix *S) {
  int inc = 1;
  int N = A->N;
  int Nm1 = N-1;
  //S->N=Nm1;
  //assert(S->N==N-1);
  assert(Nm1>0);
  copySubMatrix(A,S,Nm1);

  double factor = -1./ELEM(A,Nm1,Nm1);
  // this next line does S = S - A12 A22^(-1) A21;
  dger_(&S->N, &S->N, &factor, &ELEM(A,0,Nm1), &inc, &ELEM(A,Nm1,0), &A->N, S->data, &S->N);
  return 0;
}



int addRowColToInverse(matrix const*A, vector const*Rtilde, vector const*Qtilde, double const sTilde, matrix *Ap1) {
  int inc = 1;
  int N = A->N;
  assert(N==Rtilde->N);
  assert(N==Qtilde->N);
  int Np1 = N+1;
  copySubMatrix(A,Ap1,Np1);

  double one  = 1.0;
  double minus=-1.0;
  double zero = 0.0;
  
  printf ("\nsTilde=%f\n",sTilde);
  double factor = 1.0/sTilde;
  // got pAccept now.
  dger_(&N, &N, &factor, Qtilde->data, &inc, Rtilde->data, &inc, Ap1->data, &Np1);
  printf ("\nAp1=\n"); printMatrix(Ap1);
  
  //daxpy_(&N, &factor, Qtilde->data, &inc, &ELEM(Ap1,0,N), &inc);				
  //printf ("\nAp1=\n"); printMatrix(Ap1);
  
  dcopy_(&N, Rtilde->data, &inc, &ELEM(Ap1,N,0), &Np1);
  dcopy_(&N, Qtilde->data, &inc, &ELEM(Ap1,0,N), &inc);
  
  
  printf ("\nAp1=\n"); printMatrix(Ap1);
  
  ELEM(Ap1,N,N)=sTilde;
  
  
  
  

  //double factor = -1./ELEM(A,Nm1,Nm1);
  // this next line does S = S - A12 A22^(-1) A21;
  //dger_(&S->N, &S->N, &factor, &ELEM(A,0,Nm1), &inc, &ELEM(A,Nm1,0), &A->N, S->data, &S->N);
  return 0;
}
		































// ---------------------------- testing routines -------------------------------



// adding the vector U and V to the matrix A plus the corner should, in principle
// give the same matrix as if you invert, apply shermanMorrison and invert again.
int testAddRowColToInverse(int verbose) {
  matrix A,B,Sol;
  resizeMatrix(&A,3);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  transposeMatrix(&A);
  
  resizeMatrix(&Sol,4);
  p = Sol.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0; *p++ = 5.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0; *p++ = 1.0;
  *p++ = 1.0; *p++ = 3.0; *p++ = 5.0; *p++ = 6.0;
  transposeMatrix(&Sol);
  
  if(verbose) {
    printf("\nA=\n"); printMatrix(&A); 
    printf("\nSol=\n"); printMatrix(&Sol);
  }
  invertMatrix(&Sol);
  if(verbose) {printf("\ninverting\nSol=\n"); printMatrix(&Sol);}
  
  vector Q, R, Qtilde, Rtilde;
  resizeVector(&Q,3);
  resizeVector(&R,3);
  resizeVector(&Qtilde,3);
  resizeVector(&Rtilde,3);
  p = Q.data; *p++ = 5.0; *p++ = 2.0; *p++ = 1.0;
  p = R.data; *p++ = 1.0; *p++ = 3.0; *p++ = 5.0;
  
  double s = 6.0;
  
  invertMatrix(&A);
  if(verbose) {printf("\ninverting\nA=\n"); printMatrix(&A);}
  if(verbose) {
    printf("\nQ=\n"); printVector(&Q); 
    printf("\nR=\n"); printVector(&R);
  }
  matrixVectorProduct(&A, &Q, 1.0,&Qtilde);
  double sTilde = 1.0/( s - scalarProduct(&R,&Qtilde) );
  vectorMatrixProduct(&R, &A, -sTilde,&Rtilde);
  scaleVector(&Qtilde,-sTilde);
  
  printf("sTilde=%f",sTilde);
  if(verbose) {
    printf("\nQtilde=\n"); printVector(&Qtilde); 
    printf("\nRtilde=\n"); printVector(&Rtilde);
  }
  
  addRowColToInverse(&A,&Rtilde,&Qtilde,sTilde,&B);
  if(verbose) {
    printf("\ncalculating\nB=\n"); printMatrix(&B);
    printf("%s\n",areEqual_M(&Sol,&B)? "B==Sol": "B!=Sol");
  }
  
  return !areEqual_M(&Sol,&B);
}


int testSchurComplement(int verbose) {
  matrix A,B,S;
  resizeMatrix(&A,4);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0; *p++ = 0.9;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0; *p++ = 0.3;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0; *p++ = 8.9;
  *p++ = 2.0; *p++ = 1.5; *p++ = 5.0; *p++ = 3.0;
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  
  schurComplement(&A,&S);
  if(verbose) {printf("\ncalculate Schur complement\nS=\n"); printMatrix(&S);}

  invertMatrix(&A);
  if(verbose) {printf("\nA^-1=\n"); printMatrix(&A);}
  copySubMatrix(&A,&B,3);

  invertMatrix(&B);
  if(verbose) {
    printf("\ncalculate inverse of the 3x3 submatrix of A^-1\nB=\n"); printMatrix(&B);
    printf("%s\n",areEqual_M(&S,&B)? "B==S": "B!=S");
  }
  
  return !areEqual_M(&S,&B);
}



int testMatrixVectorProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntestMatrixVectorProduct():\n");
  matrix A;
  resizeMatrix(&A,3);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  
  vector X;
  resizeVector(&X,3);
  vector Y;
  resizeVector(&Y,3);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  if(verbose) {printf("\nX=\n"); printVector(&X);}
  
  matrixVectorProduct(&A,&X, 1.0,&Y);
  if(verbose) {printf("\nmultiplication\nY=A*X=\n"); printVector(&Y);}
  
  
  vector Sol;
  resizeVector(&Sol,3);
  p = Sol.data;
  *p++ = 29.2; *p++ = 18.3; *p++ = 13.6;
  return !areEqual_V(&Y,&Sol);
}


int testTranspose(int verbose) {
  //if(verbose) printf("\n-------------\ntestTranspose():\n");
  matrix A, Sol;
  resizeMatrix(&A,3);
  resizeMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  
  transposeMatrix(&A);
  if(verbose) {printf("\ntransposing\nA=\n"); printMatrix(&A);}
  
  p = Sol.data;
  *p++ = 1.0; *p++ = 4.0; *p++ = 3.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  return !areEqual_M(&Sol,&A);
}

int testCopy(int verbose) {
  //if(verbose) printf("\n-------------\ntestCopy():");
  matrix A, B, C;
  resizeMatrix(&A,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 4.0; *p++ = 0.0;
  *p++ = 2.0; *p++ = 1.5; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  transposeMatrix(&A); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); printMatrix(&A);
    printf("\nB=\n"); printMatrix(&B);
    printf("%s\n",areEqual_M(&A,&B)? "B==A": "B!=A");
  }

  copyMatrix(&A,&B);
  
  if(verbose){
    printf("\ncopying\nB=\n"); printMatrix(&B);
    printf("%s\n",areEqual_M(&A,&B)? "B==A": "B!=A");
  }

  copySubMatrix(&A,&C,2);
  
  if(verbose){
    printf("\ncopying\nC=\n"); printMatrix(&C);
  }

  // ----------------------------

  vector X,Y;
  resizeVector(&X,3);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  
  if(verbose){
    printf("\nX=\n"); printVector(&X);
    printf("\nY=\n"); printVector(&Y);
    printf("%s\n",areEqual_V(&X,&Y)? "X==Y": "X!=Y");
  }

  copyVector(&X,&Y);
  if(verbose){
    printf("\ncopying\nY=\n"); printVector(&Y);
    printf("%s\n",areEqual_V(&X,&Y)? "X==Y": "X!=Y");
  }

  
  return !areEqual_V(&X,&Y) + !areEqual_M(&A,&B);
}

int testMultiply(int verbose) {
  //if(verbose) printf("\n-------------\ntestMultiply():\n");
  matrix A, B, C, Sol;
  
  resizeMatrix(&A,3);
  resizeMatrix(&B,3);
  resizeMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  transposeMatrix(&A); // with this meth of input, we need to transpose
  
  p = B.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 2.0; *p++ = 2.0;
  transposeMatrix(&B); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); printMatrix(&A);
    printf("\nB=\n"); printMatrix(&B);
  }

  matrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {
    printf("\nC=A*B=\n"); printMatrix(&C);
  }

  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 4.0; *p++ = 13.0;
  transposeMatrix(&Sol); // with this meth of input, we need to transpose
  return !areEqual_M(&C,&Sol);
}

int testInvert(int verbose) {
  //if(verbose) printf("\n-------------\ntestInvert():\n");
  matrix A, B, C, Sol;
  
  resizeMatrix(&A,3);
  resizeMatrix(&B,3);
  resizeMatrix(&C,3);
  resizeMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  *p++ = 1.0; *p++ = 1.0; *p++ = 4.0;
  *p++ = 2.0; *p++ = 2.0; *p++ = 1.0;
  transposeMatrix(&A); // with this meth of input, we need to transpose
  
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  copyMatrix(&A,&B);
  if(verbose) {printf("\ncopying\nB=\n"); printMatrix(&B);}

  invertMatrix(&B);
  if(verbose) {printf("\ninverting\nB=\n"); printMatrix(&B);}
  
  matrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {printf("\nC=A*B=\n"); printMatrix(&C);}
  
  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  transposeMatrix(&Sol); // with this meth of input, we need to transpose
  return !areEqual_M(&C,&Sol);
}

int testSwaps(int verbose) {
  //if(verbose) printf("\n-------------\ntestSwapRow():\n");
  matrix A, Sol;
  resizeMatrix(&A,4);
  resizeMatrix(&Sol,4);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0; *p++ = 3.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0; *p++ = 3.0;
  transposeMatrix(&A); // with this meth of input, we need to transpose
  
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  matrixSwapRows(&A,1,3);
  if(verbose) {printf("\nswaping rows 1-3\nA=\n"); printMatrix(&A);}
  matrixSwapCols(&A,0,2);
  if(verbose) {printf("\nswaping cols 0-2\nA=\n"); printMatrix(&A);}
  
  //solution:
  p = Sol.data;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 1.0; *p++ = 2.0; *p++ = 0.0; *p++ = 3.0;
  *p++ = 1.0; *p++ = 2.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0; *p++ = 3.0;
  transposeMatrix(&Sol); // with this meth of input, we need to transpose
  return !areEqual_M(&A,&Sol);
}


int testScalarProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntestCopy():");

  vector X,Y;
  resizeVector(&X,3);
  resizeVector(&Y,3);
  double * p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  p = Y.data;
  *p++ = 2.0; *p++ = 1.3; *p++ = 0.3;
  
  if(verbose){
    printf("\nX=\n"); printVector(&X);
    printf("\nY=\n"); printVector(&Y);
  }
  
  double a=scalarProduct(&X,&Y);
  if(verbose) printf("\nscalar product\na=X.Y=% 4.4f\n", a);
  
  return !doubleEqual(a,7.79);
}


// ---------------------------- main -------------------------------------------

#define COLORNORMAL "\x1B[0m"
#define COLORRED    "\033[1m\x1B[31m"
#define COLORGREEN  "\033[1m\x1B[32m"

int passOrFail(char * fctName, int numErr)
{ 
  int numChar = strlen(fctName), N=50, i;
  printf("%s()  ",fctName);
  for(i=0;i<N-numChar;i++) printf("-");
  if(numErr==0) {
    printf("   %sPASS%s\n",COLORGREEN,COLORNORMAL); //for the color
    return 0;
  }
  else {
    printf("   %sFAIL  %d errors. %s\n",COLORRED,numErr,COLORNORMAL); 
    return 1;
  }
}

int main() {
  int verbose=0;
  int Nfail=0;
  Nfail+= passOrFail("testCopy",                testCopy(verbose));
  Nfail+= passOrFail("testMultiply",            testMultiply(verbose));
  Nfail+= passOrFail("testSwaps",               testSwaps(verbose));
  Nfail+= passOrFail("testTranspose",           testTranspose(verbose));
  Nfail+= passOrFail("testInvert",              testInvert(verbose));
  Nfail+= passOrFail("testMatrixVectorProduct", testMatrixVectorProduct(verbose));
  Nfail+= passOrFail("testScalarProduct",       testScalarProduct(verbose));
  Nfail+= passOrFail("testSchurComplement",     testSchurComplement(verbose));
  Nfail+= passOrFail("testAddRowColToInverse",  testAddRowColToInverse(1));

  if(Nfail==0) printf("%s100%% of the test PASSED%s\n\n",COLORGREEN,COLORNORMAL);
  else printf("%sOh no. ABORT!%s\n\n",COLORRED,COLORNORMAL);
  
  return 0;
}
