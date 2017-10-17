
#include "../test.h"

#include "array_test.h"

int main() {
  int verbose=1;
  int Nfail=0;
  
  printf("\n\ntesting dynamic arrays:\n\n");
  Nfail+= passOrFail("test_Array_double",       test_Array_double(verbose));
  Nfail+= passOrFail("test_Array_int",          test_Array_int(verbose));

  
  
  if(Nfail==0) printf("%s100%% of the test PASSED%s\n\n",COLORGREEN,COLORNORMAL);
  else printf("%sOh no. ABORT!%s\n\n",COLORRED,COLORNORMAL);
  
  return 0;
}
