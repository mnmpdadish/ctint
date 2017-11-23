

#include "../util/test.h"

#include "../util/stringUtil_test.h"
#include "dMatrix_test.h"
#include "cMatrix_test.h"
#include "../findGreenSymmetries/findGreenSymmetries_test.h"
#include "../oneBodyMatrix/oneBodyMatrix_test.h"

#include "../util/arrays/array_test.h"
#include "../util/arrays/array.h"

int main() {
  int verbose=1;
  int Nfail=0;
  
  printf("\n\n--------------------\ntesting matrix double operations:\n--------------------\n");
  Nfail+= passOrFail("test_dCopy",                  test_dCopy(verbose));
  Nfail+= passOrFail("test_dMultiply",              test_dMultiply(verbose));
  Nfail+= passOrFail("test_dSwaps",                 test_dSwaps(verbose));
  Nfail+= passOrFail("test_dTranspose",             test_dTranspose(verbose));
  Nfail+= passOrFail("test_dInvert",                test_dInvert(verbose));
  Nfail+= passOrFail("test_dMatrixVectorProduct",   test_dMatrixVectorProduct(verbose));
  Nfail+= passOrFail("test_dScalarProduct",         test_dScalarProduct(verbose));
  Nfail+= passOrFail("test_dSchurComplement",       test_dSchurComplement(verbose));        // used for remove vertex
  Nfail+= passOrFail("test_dAddRowColToInverse",    test_dAddRowColToInverse(verbose));     // used for insert vertex
  Nfail+= passOrFail("test_dAddOneElementToInverse",test_dAddOneElementToInverse(verbose));  // used for spin-flip
  

  
  printf("\n\n--------------------\ntesting matrix complex operations:\n--------------------\n");
  Nfail+= passOrFail("test_cCopy",                test_cCopy(verbose));
  Nfail+= passOrFail("test_cMultiply",            test_cMultiply(verbose));
  Nfail+= passOrFail("test_cDag",                 test_cDag(verbose));
  Nfail+= passOrFail("test_cInvert",              test_cInvert(verbose));
  Nfail+= passOrFail("test_cMatrixVectorProduct", test_cMatrixVectorProduct(verbose));
  Nfail+= passOrFail("test_cAddition",            test_cAddition(verbose));

  //Nfail+= passOrFail("test_huge_dMatrix",            test_huge_dMatrix(verbose));

  verdict(Nfail);
  
  return 0;
}
