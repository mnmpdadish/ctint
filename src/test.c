

#include "util/test.h"

#include "util/stringUtil_test.h"
#include "matrix/dMatrix_test.h"
#include "matrix/cMatrix_test.h"
#include "findGreenSymmetries/findGreenSymmetries_test.h"
#include "oneBodyMatrix/oneBodyMatrix_test.h"

#include "util/arrays/array_test.h"
#include "util/arrays/array.h"

int main() {
  int verbose=0;
  int Nfail=0;
  printf("\n\n--------------------\ntesting dynamic arrays:\n--------------------\n");
  Nfail+= passOrFail("test_Array_double",       test_Array_double(verbose));
  Nfail+= passOrFail("test_Array_int",          test_Array_int(verbose));

  printf("\n\n\n--------------------\ntesting string utils:\n--------------------\n");
  Nfail+= passOrFail("test_readIntInOneLine",   test_readIntInOneLine(verbose));

  printf("\n\n--------------------\ntesting matrix double operations:\n--------------------\n");
  Nfail+= passOrFail("test_dCopy",                test_dCopy(verbose));
  Nfail+= passOrFail("test_dMultiply",            test_dMultiply(verbose));
  Nfail+= passOrFail("test_dSwaps",               test_dSwaps(verbose));
  Nfail+= passOrFail("test_dTranspose",           test_dTranspose(verbose));
  Nfail+= passOrFail("test_dInvert",              test_dInvert(verbose));
  Nfail+= passOrFail("test_dMatrixVectorProduct", test_dMatrixVectorProduct(verbose));
  Nfail+= passOrFail("test_dScalarProduct",       test_dScalarProduct(verbose));
  Nfail+= passOrFail("test_dSchurComplement",     test_dSchurComplement(verbose));
  Nfail+= passOrFail("test_dAddRowColToInverse",  test_dAddRowColToInverse(verbose));

  
  printf("\n\n--------------------\ntesting matrix complex operations:\n--------------------\n");
  Nfail+= passOrFail("test_cCopy",                test_cCopy(verbose));
  Nfail+= passOrFail("test_cMultiply",            test_cMultiply(verbose));
  Nfail+= passOrFail("test_cDag",                 test_cDag(verbose));
  Nfail+= passOrFail("test_cInvert",              test_cInvert(verbose));
  Nfail+= passOrFail("test_cMatrixVectorProduct", test_cMatrixVectorProduct(verbose));
  Nfail+= passOrFail("test_cAddition",            test_cAddition(verbose));

  printf("\n\n--------------------\ntesting Green function \n--------------------\n");
  Nfail+= passOrFail("test_readSymmetriesPlaquettes2x2",       test_readSymmetriesPlaquettes2x2(verbose,"files/plaquette2x2.model"));
  Nfail+= passOrFail("test_readSymmetriesPlaquettes4x4",       test_readSymmetriesPlaquettes4x4(verbose,"files/plaquette4x4.model"));

  printf("\n\n--------------------\ntesting one body matrix \n--------------------\n");
  Nfail+= passOrFail("test_oneBodyMatrix",       test_oneBodyMatrix(verbose,"files/plaquette2x2.model"));
  

  verdict(Nfail);
  
  return 0;
}
