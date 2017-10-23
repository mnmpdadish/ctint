
#include "../util/test.h"
#include "oneBodyMatrix_test.h"

int main() {
  int verbose=1;
  int Nfail=0;
  printf("\n--------------------\ntesting one body matrix \n--------------------\n\n");
  Nfail+= passOrFail("test_oneBodyMatrix",       test_oneBodyMatrix(verbose));
  
  verdict(Nfail);
  
  return 0;
}

