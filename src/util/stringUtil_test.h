#include "stringUtil.h"



int test_readIntInOneLine(int verbose)
{ 
  //test 1
  char line[] = "1 3 2 \t \t\t  34784 5 #this is a comment 4 23";
  unsigned int nElement = countElementInStr(line, " \t\n");
  int *arrayInt = (int *) malloc(nElement * sizeof(int)); 
  readIntInStr(line, arrayInt, " \t\n");
  int i, Nerror=0;
  int solution1[5]={1,3,2,34784,5};
  
  if(verbose){
    printf("integer read:\n");
    for(i=0;i<nElement;i++) printf("%d ", arrayInt[i]);
    printf("\n");
  }
  
  for(i=0;i<nElement;i++) if(solution1[i] != arrayInt[i]) Nerror++;
  free(arrayInt);
  
  //test 2
  char parenthesisTuple[] = "(1,2,4)";
  int tupleInt[3] = {0,0,0};
  readIntInStr(parenthesisTuple, tupleInt, "(,)\t\n");
  int solution2[3]={1,2,4};
  
  if(verbose){
    printf("integer read:\n");
    for(i=0;i<3;i++) printf("%d ", tupleInt[i]);
    printf("\n");
  }
  
  for(i=0;i<3;i++) if(solution2[i] != tupleInt[i]) Nerror++;
  
  return Nerror;
}


/*
int test_readIntInParenthesis(int verbose)
{
  char parenthesisTuple[] = "(1,2,4)";
  int arrayInt[3] = {0,0,0};
  readIntInStr(parenthesisTuple, arrayInt, "(,)\t\n");
  //readIntInParenthesis(parenthesisTuple,arrayInt);
  int i, Nerror=0;
  int solution[5]={1,2,4};
  
  if(verbose){
    printf("integer read:\n");
    for(i=0;i<3;i++) printf("%d ", arrayInt[i]);
    printf("\n");
  }
  
  for(i=0;i<3;i++) if(solution[i] != arrayInt[i]) Nerror++;
  
  return Nerror;
}*/


