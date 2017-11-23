#pragma once

#include "dMatrix.h"

// ---------------------------- testing routines -------------------------------




int test_dCopy(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dCopy():");
  dMatrix A, B, C;
  init_dMatrix(&A,3);
  init_dMatrix(&B,0);
  init_dMatrix(&C,0);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 4.0; *p++ = 0.0;
  *p++ = 2.0; *p++ = 1.5; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  transpose_dMatrix(&A); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); print_dMatrix(&A);
    printf("\nB=\n"); print_dMatrix(&B);
    printf("%s\n",areEqual_dMatrix(&A,&B)? "B==A": "B!=A");
  }

  copy_dMatrix(&A,&B);
  
  if(verbose){
    printf("\ncopying\nB=\n"); print_dMatrix(&B);
    printf("%s\n",areEqual_dMatrix(&A,&B)? "B==A": "B!=A");
  }

  copySub_dMatrix(&A,&C,2);
  
  if(verbose){
    printf("\ncopying\nC=\n"); print_dMatrix(&C);
  }

  // ----------------------------

  dVector X,Y;
  init_dVector(&X,3);
  init_dVector(&Y,0);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  
  if(verbose){
    printf("\nX=\n"); print_dVector(&X);
    printf("\nY=\n"); print_dVector(&Y);
    printf("%s\n",areEqual_dVector(&X,&Y)? "X==Y": "X!=Y");
  }

  copy_dVector(&X,&Y);
  if(verbose){
    printf("\ncopying\nY=\n"); print_dVector(&Y);
    printf("%s\n",areEqual_dVector(&X,&Y)? "X==Y": "X!=Y");
  }

  int Nerror = !areEqual_dVector(&X,&Y) + !areEqual_dMatrix(&A,&B);
  free_dMatrix(&A);
  free_dMatrix(&B);
  free_dMatrix(&C);
  free_dVector(&X);
  free_dVector(&Y);
  return Nerror;
}



// adding the Vector U and V to the matrix A plus the corner should, in principle
// give the same matrix as if you invert, apply "shermanMorrison" and invert again.
int test_dAddRowColToInverse(int verbose) {
  dMatrix A,B,Sol;
  init_dMatrix(&A,3);
  init_dMatrix(&B,3);
  init_dMatrix(&Sol,4);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  transpose_dMatrix(&A);
  
  p = Sol.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0; *p++ = 5.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0; *p++ = 1.0;
  *p++ = 1.0; *p++ = 3.0; *p++ = 5.0; *p++ = 6.0;
  transpose_dMatrix(&Sol);
  
  if(verbose) {
    printf("\nA=\n"); print_dMatrix(&A); 
    printf("\nSol=\n"); print_dMatrix(&Sol);
  }
  invert_dMatrix(&Sol);
  if(verbose) {printf("\ninverting\nSol=\n"); print_dMatrix(&Sol);}
  
  dVector Q, R, Qtilde, Rtilde;
  init_dVector(&Q,3);
  init_dVector(&R,3);
  init_dVector(&Qtilde,3);
  init_dVector(&Rtilde,3);
  p = Q.data; *p++ = 5.0; *p++ = 2.0; *p++ = 1.0;
  p = R.data; *p++ = 1.0; *p++ = 3.0; *p++ = 5.0;
  
  double s = 6.0;
  
  invert_dMatrix(&A);
  if(verbose) {printf("\ninverting\nA=\n"); print_dMatrix(&A);}
  if(verbose) {
    printf("\nQ=\n"); print_dVector(&Q); 
    printf("\nR=\n"); print_dVector(&R);
  }
  dMatrixVectorProduct(&A, &Q, 1.0,&Qtilde);
  double sTilde = 1.0/( s - dScalarProduct(&R,&Qtilde) );
  dVectorMatrixProduct(&R, &A, -sTilde,&Rtilde);
  scale_dVector(&Qtilde,-sTilde);
  
  //printf("sTilde=%f",sTilde);
  if(verbose) {
    printf("\nQtilde=\n"); print_dVector(&Qtilde); 
    printf("\nRtilde=\n"); print_dVector(&Rtilde);
  }
  
  dAddRowColToInverse(&A,&Rtilde,&Qtilde,sTilde,&B);
  if(verbose) {
    printf("\ncalculating\nB=\n"); print_dMatrix(&B);
    printf("%s\n",areEqual_dMatrix(&Sol,&B)? "B==Sol": "B!=Sol");
  }
  

  int Nerror = !areEqual_dMatrix(&Sol,&B);
  free_dMatrix(&A);
  free_dMatrix(&B);
  free_dVector(&Q);
  free_dVector(&R);
  free_dVector(&Qtilde);
  free_dVector(&Rtilde);
  return Nerror;
}


int test_dSchurComplement(int verbose) {
  dMatrix A,B,S;
  init_dMatrix(&A,4);
  init_dMatrix(&B,4);
  init_dMatrix(&S,4);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0; *p++ = 0.9;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0; *p++ = 0.3;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0; *p++ = 8.9;
  *p++ = 2.0; *p++ = 1.5; *p++ = 5.0; *p++ = 3.0;
  if(verbose) {printf("\nA=\n"); print_dMatrix(&A);}
  
  dSchurComplement(&A,&S);
  if(verbose) {printf("\ncalculate Schur complement\nS=\n"); print_dMatrix(&S);}

  invert_dMatrix(&A);
  if(verbose) {printf("\nA^-1=\n"); print_dMatrix(&A);}
  copySub_dMatrix(&A,&B,3);

  invert_dMatrix(&B);
  if(verbose) {
    printf("\ncalculate inverse of the 3x3 submatrix of A^-1\nB=\n"); print_dMatrix(&B);
    printf("%s\n",areEqual_dMatrix(&S,&B)? "B==S": "B!=S");
  }
  
  int Nerror = !areEqual_dMatrix(&S,&B);
  free_dMatrix(&A);
  free_dMatrix(&B);
  free_dMatrix(&S);
  return Nerror;
}



int test_dMatrixVectorProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dMatrixVectorProduct():\n");
  dMatrix A;
  init_dMatrix(&A,3);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); print_dMatrix(&A);}
  
  dVector X;
  init_dVector(&X,3);
  dVector Y;
  init_dVector(&Y,3);
  p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  if(verbose) {printf("\nX=\n"); print_dVector(&X);}
  
  dMatrixVectorProduct(&A,&X, 1.0,&Y);
  if(verbose) {printf("\nmultiplication\nY=A*X=\n"); print_dVector(&Y);}
  
  
  dVector Sol;
  init_dVector(&Sol,3);
  p = Sol.data;
  *p++ = 29.2; *p++ = 18.3; *p++ = 13.6;

  int Nerror = !areEqual_dVector(&Y,&Sol);
  free_dMatrix(&A);
  free_dVector(&X);
  free_dVector(&Y);
  free_dVector(&Sol);
  return Nerror;
}


int test_dTranspose(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dTranspose():\n");
  dMatrix A, Sol;
  init_dMatrix(&A,3);
  init_dMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  if(verbose) {printf("\nA=\n"); print_dMatrix(&A);}
  
  transpose_dMatrix(&A);
  if(verbose) {printf("\ntransposing\nA=\n"); print_dMatrix(&A);}
  
  p = Sol.data;
  *p++ = 1.0; *p++ = 4.0; *p++ = 3.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  
  int Nerror = !areEqual_dMatrix(&Sol,&A);
  free_dMatrix(&A);
  free_dMatrix(&Sol);
  return Nerror;
}


int test_dMultiply(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dMultiply():\n");
  dMatrix A, B, C, Sol;
  
  init_dMatrix(&A,3);
  init_dMatrix(&B,3);
  init_dMatrix(&C,0);
  init_dMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0;
  transpose_dMatrix(&A); // with this meth of input, we need to transpose
  
  p = B.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 2.0; *p++ = 2.0;
  transpose_dMatrix(&B); // with this meth of input, we need to transpose

  if(verbose){
    printf("\nA=\n"); print_dMatrix(&A);
    printf("\nB=\n"); print_dMatrix(&B);
  }

  dMatrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {
    printf("\nC=A*B=\n"); print_dMatrix(&C);
  }

  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 5.5;
  *p++ = 0.0; *p++ = 4.0; *p++ = 13.0;
  transpose_dMatrix(&Sol); // with this meth of input, we need to transpose
  
  int Nerror = !areEqual_dMatrix(&C,&Sol);
  free_dMatrix(&A);
  free_dMatrix(&B);
  free_dMatrix(&C);
  free_dMatrix(&Sol);
  return Nerror;
}

int test_dInvert(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dInvert():\n");
  dMatrix A, B, C, Sol;
  
  init_dMatrix(&A,3);
  init_dMatrix(&B,3);
  init_dMatrix(&C,3);
  init_dMatrix(&Sol,3);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  *p++ = 1.0; *p++ = 1.0; *p++ = 4.0;
  *p++ = 2.0; *p++ = 2.0; *p++ = 1.0;
  transpose_dMatrix(&A); // with this meth of input, we need to transpose
  
  if(verbose) {printf("\nA=\n"); print_dMatrix(&A);}
  copy_dMatrix(&A,&B);
  if(verbose) {printf("\ncopying\nB=\n"); print_dMatrix(&B);}

  invert_dMatrix(&B);
  if(verbose) {printf("\ninverting\nB=\n"); print_dMatrix(&B);}
  
  dMatrixMatrixMultiplication(&A,&B,&C);
  if(verbose) {printf("\nC=A*B=\n"); print_dMatrix(&C);}
  
  //solution:
  p = Sol.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0;
  transpose_dMatrix(&Sol); // with this meth of input, we need to transpose
  
  int Nerror = !areEqual_dMatrix(&C,&Sol);
  free_dMatrix(&A);
  free_dMatrix(&B);
  free_dMatrix(&C);
  free_dMatrix(&Sol);
  return Nerror;
}

int test_dSwaps(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dSwapRow():\n");
  dMatrix A, Sol;
  init_dMatrix(&A,4);
  init_dMatrix(&Sol,4);
  
  double * p = A.data;
  *p++ = 1.0; *p++ = 0.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0; *p++ = 3.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 2.0; *p++ = 1.0; *p++ = 3.0;
  transpose_dMatrix(&A); // with this meth of input, we need to transpose
  
  if(verbose) {printf("\nA=\n"); print_dMatrix(&A);}
  dMatrixSwapRows(&A,1,3);
  if(verbose) {printf("\nswaping rows 1-3\nA=\n"); print_dMatrix(&A);}
  dMatrixSwapCols(&A,0,2);
  if(verbose) {printf("\nswaping cols 0-2\nA=\n"); print_dMatrix(&A);}
  
  //solution:
  p = Sol.data;
  *p++ = 0.0; *p++ = 0.0; *p++ = 1.0; *p++ = 0.0;
  *p++ = 1.0; *p++ = 2.0; *p++ = 0.0; *p++ = 3.0;
  *p++ = 1.0; *p++ = 2.0; *p++ = 0.0; *p++ = 0.0;
  *p++ = 0.0; *p++ = 1.0; *p++ = 0.0; *p++ = 3.0;
  transpose_dMatrix(&Sol); // with this meth of input, we need to transpose
  
  int Nerror = !areEqual_dMatrix(&A,&Sol);
  free_dMatrix(&A);
  free_dMatrix(&Sol);
  return Nerror;
}


int test_dScalarProduct(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dCopy():");

  dVector X,Y;
  init_dVector(&X,3);
  init_dVector(&Y,3);
  double * p = X.data;
  *p++ = 1.0; *p++ = 3.3; *p++ = 5.0;
  p = Y.data;
  *p++ = 2.0; *p++ = 1.3; *p++ = 0.3;
  
  if(verbose){
    printf("\nX=\n"); print_dVector(&X);
    printf("\nY=\n"); print_dVector(&Y);
  }
  
  double a=dScalarProduct(&X,&Y);
  if(verbose) printf("\nscalar product\na=X.Y=% 4.4f\n", a);
    
  int Nerror = !doubleEqual(a,7.79);
  free_dVector(&X);
  free_dVector(&Y);
  return Nerror;
}




/*
int test_huge_dMatrix(int verbose) {
  //if(verbose) printf("\n-------------\ntest_dCopy():");
  dMatrix A;
  init_dMatrix(&A,12018);
  //init_dMatrix(&B,0);
  free_dMatrix(&A);
  //free_dMatrix(&B);
  return 0;
}
*/












// adding the Vector U and V to the matrix A plus the corner should, in principle
// give the same matrix as if you invert, apply "shermanMorrison" and invert again.
int test_dAddOneElementToInverse(int verbose) {
  dMatrix A, Sol;
  init_dMatrix(&A,3);
  init_dMatrix(&Sol,3);
  double * p = A.data;
  *p++ = 1.0; *p++ = 5.0; *p++ = 2.0;
  *p++ = 4.0; *p++ = 1.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0; *p++ = 1.0;
  transpose_dMatrix(&A);
  
  p = Sol.data;
  *p++ = 1.0; *p++ = 5.0;  *p++ = 2.0;
  *p++ = 4.0; *p++ =-29.0; *p++ = 2.0;
  *p++ = 3.0; *p++ = 2.0;  *p++ = 1.0;
  transpose_dMatrix(&Sol);
  
  if(verbose) {
    printf("\nA=\n"); print_dMatrix(&A); 
    printf("\nSol=\n"); print_dMatrix(&Sol);
  }
  
  invert_dMatrix(&A);
  if(verbose) {
    printf("\ninverting\nA=\n"); 
    print_dMatrix(&A);
  }
  
  dVector Row, Col;
  init_dVector(&Row,0);
  init_dVector(&Col,0);
  dCopyRowIntoVector(&A, &Row, 1);
  dCopyColIntoVector(&A, &Col, 1);

  if(verbose) {
    printf("\nRow=\n"); print_dVector(&Row); 
    printf("\nCol=\n"); print_dVector(&Col);
  }
  
  double value=-30.;
  double factor=-value/(1+value*ELEM_VAL(A,1,1));
  if(verbose) printf("value=%f,  factor=%f\n",value, factor);
  dAddOneElementToInverse(&A, &Row, &Col, factor);
  
  invert_dMatrix(&A);
  if(verbose) {printf("\ninverting\nA=\n"); print_dMatrix(&A);}
  
  int Nerror = !areEqual_dMatrix(&Sol,&A);
  free_dMatrix(&A);
  free_dMatrix(&Sol);
  free_dVector(&Row);
  free_dVector(&Col);
  return Nerror;
}

