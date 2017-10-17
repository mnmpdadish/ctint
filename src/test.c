
#include "util/test.h"

#include "findGreenTools/stringUtil_test.h"
#include "matrix_test.h"

#include "util/arrays/array_test.h"
#include "util/arrays/array.h"

int main() {
  int verbose=0;
  int Nfail=0;
  printf("\n\ntesting matrix operations:\n\n");
  Nfail+= passOrFail("testCopy",                testCopy(verbose));
  Nfail+= passOrFail("testMultiply",            testMultiply(verbose));
  Nfail+= passOrFail("testSwaps",               testSwaps(verbose));
  Nfail+= passOrFail("testTranspose",           testTranspose(verbose));
  Nfail+= passOrFail("testInvert",              testInvert(verbose));
  Nfail+= passOrFail("testMatrixVectorProduct", testMatrixVectorProduct(verbose));
  Nfail+= passOrFail("testScalarProduct",       testScalarProduct(verbose));
  Nfail+= passOrFail("testSchurComplement",     testSchurComplement(verbose));
  Nfail+= passOrFail("testAddRowColToInverse",  testAddRowColToInverse(verbose));

  verbose=1;
  printf("\n\ntesting string utils:\n\n");
  Nfail+= passOrFail("test_readIntInOneLine",   test_readIntInOneLine(verbose));
  
  printf("\n\ntesting dynamic arrays:\n\n");
  Nfail+= passOrFail("test_Array_double",       test_Array_double(verbose));
  Nfail+= passOrFail("test_Array_int",          test_Array_int(verbose));

  
  
  if(Nfail==0) printf("%s100%% of the test PASSED%s\n\n",COLORGREEN,COLORNORMAL);
  else printf("%sOh no. ABORT!%s\n\n",COLORRED,COLORNORMAL);
  
  return 0;
}
