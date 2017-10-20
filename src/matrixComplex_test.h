#pragma once

#include "matrixComplex.h"

// ---------------------------- testing routines -------------------------------




int test_cCopy(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cCopy():");
  cMatrix A, B, C;
  init_cMatrix(&A,3);
  init_cMatrix(&B,0);
  init_cMatrix(&C,0);
  
  double complex * p = A.data;
  *p++ = 1.0+3.4*I; *p++ = 4.0;       *p++ = 0.0+1.4*I;
  *p++ = 2.0+0.0*I; *p++ = 1.5+1.*I;  *p++ = 0.0+3.2*I;
  *p++ = 0.0+1.1*I; *p++ = 2.0+0.1*I; *p++ = 1.0;
  transpose_cMatrix(&A); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); print_cMatrix(&A);
    printf("\nB=\n"); print_cMatrix(&B);
    printf("%s\n",areEqual_cMatrix(&A,&B)? "B==A": "B!=A");
  }

  copy_cMatrix(&A,&B);
  
  if(verbose){
    printf("\ncopying\nB=\n"); print_cMatrix(&B);
    printf("%s\n",areEqual_cMatrix(&A,&B)? "B==A": "B!=A");
  }

  copySub_cMatrix(&A,&C,2);
  
  if(verbose){
    printf("\ncopying\nC=\n"); print_cMatrix(&C);
  }

  // ----------------------------

  cVector X,Y;
  init_cVector(&X,3);
  init_cVector(&Y,0);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  
  if(verbose){
    printf("\nX=\n"); print_cVector(&X);
    printf("\nY=\n"); print_cVector(&Y);
    printf("%s\n",areEqual_cVector(&X,&Y)? "X==Y": "X!=Y");
  }

  copy_cVector(&X,&Y);
  if(verbose){
    printf("\ncopying\nY=\n"); print_cVector(&Y);
    printf("%s\n",areEqual_cVector(&X,&Y)? "X==Y": "X!=Y");
  }

  int Nerror = !areEqual_cVector(&X,&Y) + !areEqual_cMatrix(&A,&B);
  free_cMatrix(&A);
  free_cMatrix(&B);
  free_cMatrix(&C);
  free_cVector(&X);
  free_cVector(&Y);
  return Nerror;
}


int test_cMatrixVectorProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cMatrixVectorProduct():\n");
  cMatrix A;
  init_cMatrix(&A,3);
  double complex * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); print_cMatrix(&A);}
  
  cVector X;
  init_cVector(&X,3);
  cVector Y;
  init_cVector(&Y,3);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  if(verbose) {printf("\nX=\n"); print_cVector(&X);}
  
  cMatrixVectorProduct(&A,&X, 1.0,&Y);
  if(verbose) {printf("\nmultiplication\nY=A*X=\n"); print_cVector(&Y);}
  
  
  cVector Sol;
  init_cVector(&Sol,3);
  p = Sol.data;
  *p++ = 29.2; *p++ = 18.3; *p++ = 13.6;

  int Nerror = !areEqual_cVector(&Y,&Sol);
  free_cMatrix(&A);
  free_cVector(&X);
  free_cVector(&Y);
  free_cVector(&Sol);
  return Nerror;
}


int test_cDag(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cTranspose():\n");
  cMatrix A, Sol;
  init_cMatrix(&A,3);
  init_cMatrix(&Sol,3);
  
  double complex * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); print_cMatrix(&A);}
  
  dag_cMatrix(&A);
  if(verbose) {printf("\ndagging\nA=\n"); print_cMatrix(&A);}
  
  p = Sol.data;
  *p++ = 1.0; *p++ = 4.0; *p++ = 3.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  
  int Nerror = !areEqual_cMatrix(&Sol,&A);
  free_cMatrix(&A);
  free_cMatrix(&Sol);
  return Nerror;
}


int test_cMultiply(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cMultiply():\n");
  cMatrix A, B, C, Sol;
  
  init_cMatrix(&A,3);
  init_cMatrix(&B,3);
  init_cMatrix(&C,0);
  init_cMatrix(&Sol,3);
  
  double complex * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  transpose_cMatrix(&A); // with this meth of input, we need to transpose
  
  p = B.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 2.0; *p++ = 2.0;
  transpose_cMatrix(&B); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); print_cMatrix(&A);
    printf("\nB=\n"); print_cMatrix(&B);
  }

  cMatrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {
    printf("\nC=A*B=\n"); print_cMatrix(&C);
  }

  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 4.0; *p++ = 13.0;
  transpose_cMatrix(&Sol); // with this meth of input, we need to transpose
  
  int Nerror = !areEqual_cMatrix(&C,&Sol);
  free_cMatrix(&A);
  free_cMatrix(&B);
  free_cMatrix(&C);
  free_cMatrix(&Sol);
  return Nerror;
}

int test_cInvert(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cInvert():\n");
  cMatrix A, B, C, Sol;
  
  init_cMatrix(&A,3);
  init_cMatrix(&B,3);
  init_cMatrix(&C,3);
  init_cMatrix(&Sol,3);
  
  double complex * p = A.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  *p++ = 1.0; *p++ = 1.0; *p++ = 4.0;
  *p++ = 2.0; *p++ = 2.0; *p++ = 1.0;
  transpose_cMatrix(&A); // with this meth of input, we need to transpose
  
  if(verbose) {printf("\nA=\n"); print_cMatrix(&A);}
  copy_cMatrix(&A,&B);
  if(verbose) {printf("\ncopying\nB=\n"); print_cMatrix(&B);}

  invert_cMatrix(&B);
  if(verbose) {printf("\ninverting\nB=\n"); print_cMatrix(&B);}
  
  cMatrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {printf("\nC=A*B=\n"); print_cMatrix(&C);}
  
  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  transpose_cMatrix(&Sol); // with this meth of input, we need to transpose
  
  int Nerror = !areEqual_cMatrix(&C,&Sol);
  free_cMatrix(&A);
  free_cMatrix(&B);
  free_cMatrix(&C);
  free_cMatrix(&Sol);
  return Nerror;
}


/*
int test_cScalarProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntest_cCopy():");

  cVector X,Y;
  init_cVector(&X,3);
  init_cVector(&Y,3);
  double complex * p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  p = Y.data;
  *p++ = 2.0; *p++ = 1.3; *p++ = 0.3;
  
  if(verbose){
    printf("\nX=\n"); print_cVector(&X);
    printf("\nY=\n"); print_cVector(&Y);
  }
  
  double complex a=cScalarProduct(&X,&Y);
  if(verbose) printf("\nscalar product\na=X.Y=% 4.4f\n", a);
    
  int Nerror = !complexEqual(a,7.79);
  free_cVector(&X);
  free_cVector(&Y);
  return Nerror;
}
*/


