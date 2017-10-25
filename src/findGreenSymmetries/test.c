
#include "../util/stringUtil.h"
#include "../util/arrays/array.h"
#include "../util/test.h"
#include "findGreenSymmetries_test.h"

int main() {
  int verbose=1;
  int Nfail=0;
  printf("\n--------------------\ntesting Green function \n--------------------\n\n");
  Nfail+= passOrFail("test_readSymmetriesPlaquettes2x2",       test_readSymmetriesPlaquettes2x2(verbose));
  Nfail+= passOrFail("test_readSymmetriesPlaquettes4x4",       test_readSymmetriesPlaquettes4x4(verbose));
  Nfail+= passOrFail("test_arbitrary",                         test_arbitrary(verbose));
  
  verdict(Nfail);
  
  return 0;
}
