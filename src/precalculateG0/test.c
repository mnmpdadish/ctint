
#include "../util/stringUtil.h"
#include "../util/arrays/array.h"
#include "../util/test.h"
#include "calculateG0_test.h"

int main() {
  int verbose=1;
  int Nfail=0;
  //printf("%f %f\n", creal(1.0/(1.+I*2.)),cimag(1.0/(1.+I*2.)));
  printf("\n--------------------\ntesting precalculate Green function \n--------------------\n\n");
  Nfail+= passOrFail("test_Plaquettes2x2",       test_Plaquettes2x2(verbose));
  
  verdict(Nfail);
  
  return 0;
}
