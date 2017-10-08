#pragma once

#include "matrix.h"

// ---------------------------- testing routines -------------------------------



// adding the vector U and V to the matrix A plus the corner should, in principle
// give the same matrix as if you invert, apply "shermanMorrison" and invert again.
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



