
#include "stringUtil_test.h"
#include "../util/test.h"

int main0()
{
  char c=0; 
  int i;
  for(i=0;i<256;i++) printf("%c",c++);
  return 0; 
} 


int main() {
  int verbose=1;
  int Nfail=0;
  Nfail+= passOrFail("test_readIntInOneLine",        test_readIntInOneLine(verbose));
  //Nfail+= passOrFail("test_readIntInParenthesis",    test_readIntInParenthesis(verbose));
  
  verdict(Nfail);
  
  return 0;
}
