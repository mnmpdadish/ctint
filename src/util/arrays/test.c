
#include "../test.h"

#include "array_test.h"

int main() {
  int verbose=1;
  int Nfail=0;
  
  printf("\n\ntesting dynamic arrays:\n\n");
  Nfail+= passOrFail("test_Array_double",       test_Array_double(verbose));
  Nfail+= passOrFail("test_Array_int",          test_Array_int(verbose));

  verdict(Nfail);
  
  return 0;
}
