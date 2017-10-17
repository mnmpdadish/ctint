#pragma once

#include "matrix.h"

// ---------------------------- testing routines -------------------------------




int testCopy(int verbose) {
  //if(verbose) printf("\n-------------\ntestCopy():");
  matrix A, B, C;
  initMatrix(&A,3);
  initMatrix(&B,0);
  initMatrix(&C,0);
  
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
  initVector(&X,3);
  initVector(&Y,0);
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

  int Nerror = !areEqual_V(&X,&Y) + !areEqual_M(&A,&B);
  freeMatrix(&A);
  freeMatrix(&B);
  freeMatrix(&C);
  freeVector(&X);
  freeVector(&Y);
  return Nerror;
}



// adding the vector U and V to the matrix A plus the corner should, in principle
// give the same matrix as if you invert, apply "shermanMorrison" and invert again.
int testAddRowColToInverse(int verbose) {
  matrix A,B,Sol;
  initMatrix(&A,3);
  initMatrix(&B,3);
  initMatrix(&Sol,4);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  transposeMatrix(&A);
  
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
  initVector(&Q,3);
  initVector(&R,3);
  initVector(&Qtilde,3);
  initVector(&Rtilde,3);
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
  
  //printf("sTilde=%f",sTilde);
  if(verbose) {
    printf("\nQtilde=\n"); printVector(&Qtilde); 
    printf("\nRtilde=\n"); printVector(&Rtilde);
  }
  
  addRowColToInverse(&A,&Rtilde,&Qtilde,sTilde,&B);
  if(verbose) {
    printf("\ncalculating\nB=\n"); printMatrix(&B);
    printf("%s\n",areEqual_M(&Sol,&B)? "B==Sol": "B!=Sol");
  }
  

  int Nerror = !areEqual_M(&Sol,&B);
  freeMatrix(&A);
  freeMatrix(&B);
  freeVector(&Q);
  freeVector(&R);
  freeVector(&Qtilde);
  freeVector(&Rtilde);
  return Nerror;
}


int testSchurComplement(int verbose) {
  matrix A,B,S;
  initMatrix(&A,4);
  initMatrix(&B,4);
  initMatrix(&S,4);
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
  
  int Nerror = !areEqual_M(&S,&B);
  freeMatrix(&A);
  freeMatrix(&B);
  freeMatrix(&S);
  return Nerror;
}



int testMatrixVectorProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntestMatrixVectorProduct():\n");
  matrix A;
  initMatrix(&A,3);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); printMatrix(&A);}
  
  vector X;
  initVector(&X,3);
  vector Y;
  initVector(&Y,3);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  if(verbose) {printf("\nX=\n"); printVector(&X);}
  
  matrixVectorProduct(&A,&X, 1.0,&Y);
  if(verbose) {printf("\nmultiplication\nY=A*X=\n"); printVector(&Y);}
  
  
  vector Sol;
  initVector(&Sol,3);
  p = Sol.data;
  *p++ = 29.2; *p++ = 18.3; *p++ = 13.6;

  int Nerror = !areEqual_V(&Y,&Sol);
  freeMatrix(&A);
  freeVector(&X);
  freeVector(&Y);
  freeVector(&Sol);
  return Nerror;
}


int testTranspose(int verbose) {
  //if(verbose) printf("\n-------------\ntestTranspose():\n");
  matrix A, Sol;
  initMatrix(&A,3);
  initMatrix(&Sol,3);
  
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
  
  int Nerror = !areEqual_M(&Sol,&A);
  freeMatrix(&A);
  freeMatrix(&Sol);
  return Nerror;
}


int testMultiply(int verbose) {
  //if(verbose) printf("\n-------------\ntestMultiply():\n");
  matrix A, B, C, Sol;
  
  initMatrix(&A,3);
  initMatrix(&B,3);
  initMatrix(&C,0);
  initMatrix(&Sol,3);
  
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
  
  int Nerror = !areEqual_M(&C,&Sol);
  freeMatrix(&A);
  freeMatrix(&B);
  freeMatrix(&C);
  freeMatrix(&Sol);
  return Nerror;
}

int testInvert(int verbose) {
  //if(verbose) printf("\n-------------\ntestInvert():\n");
  matrix A, B, C, Sol;
  
  initMatrix(&A,3);
  initMatrix(&B,3);
  initMatrix(&C,3);
  initMatrix(&Sol,3);
  
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
  
  int Nerror = !areEqual_M(&C,&Sol);
  freeMatrix(&A);
  freeMatrix(&B);
  freeMatrix(&C);
  freeMatrix(&Sol);
  return Nerror;
}

int testSwaps(int verbose) {
  //if(verbose) printf("\n-------------\ntestSwapRow():\n");
  matrix A, Sol;
  initMatrix(&A,4);
  initMatrix(&Sol,4);
  
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
  
  int Nerror = !areEqual_M(&A,&Sol);
  freeMatrix(&A);
  freeMatrix(&Sol);
  return Nerror;
}


int testScalarProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntestCopy():");

  vector X,Y;
  initVector(&X,3);
  initVector(&Y,3);
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
    
  int Nerror = !doubleEqual(a,7.79);
  freeVector(&X);
  freeVector(&Y);
  return Nerror;
}



