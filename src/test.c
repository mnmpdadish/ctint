


#include "test.h"
#include "matrix_test.h"

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
  Nfail+= passOrFail("testAddRowColToInverse",  testAddRowColToInverse(verbose));

  if(Nfail==0) printf("%s100%% of the test PASSED%s\n\n",COLORGREEN,COLORNORMAL);
  else printf("%sOh no. ABORT!%s\n\n",COLORRED,COLORNORMAL);
  
  return 0;
}
