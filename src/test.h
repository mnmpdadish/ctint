#pragma once

#include "utilities.h"


// ---------------------------- main -------------------------------------------

#define COLORNORMAL "\x1B[0m"
#define COLORRED    "\033[1m\x1B[31m"
#define COLORGREEN  "\033[1m\x1B[32m"

int passOrFail(char * fctName, int numErr)
{ 
  int numChar = strlen(fctName), N=50, i;
  printf("%s()  ",fctName);
  for(i=0;i<N-numChar;i++) printf("-");
  if(numErr==0) {
    printf("   %sPASS%s\n",COLORGREEN,COLORNORMAL); //for the color
    return 0;
  }
  else {
    printf("   %sFAIL  %d errors. %s\n",COLORRED,numErr,COLORNORMAL); 
    return 1;
  }
}
