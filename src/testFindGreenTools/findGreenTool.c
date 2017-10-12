#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>





typedef struct {
  //any cluster should not have more than 64 sites, even in 2050 xD
  unsigned int data[64];
  unsigned int nSites;
  unsigned int nAssigned;
} Permutation;


typedef struct {
  Permutation * permutation;
  int n;
} Symmetries;

int allocateSymmetries(Symmetries * sym, int nSym, int nSites) {
  sym->permutation = (Permutation *) malloc(nSym * sizeof (Permutation));
  sym->n=nSym;
  int i;
  for(i=0;i<nSym;i++){
    sym->permutation[i].nSites=nSites;
    sym->permutation[i].nAssigned=0;
  } 
  return 0;
}

int addElementToPermutation(Permutation * perm, unsigned int value) {
  unsigned int i;
  //the four next lines are there to ensure that the permutation has only one of 
  //each number, between 0 and nSites-1.
  for(i=0;i<perm->nAssigned;i++){
    assert(perm->data[i] != value);
  }
  assert(value<perm->nSites && value>=0);

  perm->data[perm->nAssigned]=value;
  perm->nAssigned++;
  assert(perm->nAssigned <= perm->nSites);
  return 0;
}


char charHexa(int digit) {
  //any cluster should not have more than 64 sites, even in 2050 xD
  char characters[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!."; 
  char char1;
  if(digit < 36 && digit >= 0) char1 = characters[digit];
  else char1 = '.';
  return char1;
}


int printSymmetries(Symmetries * sym){
  int i,j;
  for(i=0;i<sym->n;i++){
    for(j=0;j<sym->permutation[i].nAssigned;j++){
      printf("%c", charHexa(sym->permutation[i].data[j]));
    }
    printf("\n");
  } 
  return 0;
}







int countNumberOfElementInOneLine(char *tempbuff){
  char *token, *saveptr;
  const char delimiter = ' ';
  token = strtok_r(tempbuff, &delimiter, &saveptr);

  int N=0;
  while( token != NULL && *token != '\n') {
    N++;
    token = strtok_r(NULL, &delimiter, &saveptr);
  }
  return N;
}

int countLineFlag(FILE * file, char *flag) {

  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, N=0;
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        char tmpstr1[50];
        sscanf(tempbuff, "%200s\n", tmpstr1);
        if (strcmp(tmpstr1,flag)==0) { 
          found=1;
          continue;
        }
      }
      else{
        if(tempbuff[0] == '#') continue;
        else if((tempbuff[0] != '\n') && countNumberOfElementInOneLine(tempbuff) != 0) N++;
        else break;
      }
    }
    printf("\n");
  }
  return N;
}


void readSymmetries(FILE * file, int nSites, Symmetries *sym) {

  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, N, i;
  char name[]="symmetries";
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        char tmpstr1[50];
        sscanf(tempbuff, "%200s\n", tmpstr1);
        if (strcmp(tmpstr1,name)==0) { 
          found=1;
          break;
        }
      }
    }
  }
  for(i=0; i<sym->n; i++){
    if (fgets(tempbuff,256,file)) {
      char *token, *saveptr;
      const char delimiter = ' ';
      token = strtok_r(tempbuff, &delimiter, &saveptr);
      
      int N=0;
      while( token != NULL && *token != '\n') {
        N++;
        unsigned int number = atoi(token);
        addElementToPermutation(&(sym->permutation[i]),number);
        token = strtok_r(NULL, &delimiter, &saveptr);
      }
      //printf("\n %d %d",sym->permutation[i].nSites,sym->n);
      assert(sym->permutation[i].nSites == N);
    }
  }   
}


int main() {
  FILE * file = fopen("plaquette4x4.in", "rt");
  if(file == NULL) {printf("file %s not found", "plaquette4x4.in"); exit(1);}
  printf("reading symmetries from plaquette4x4.in:\n");

  int nSites=16;
  //Symmetries symmetries
  int nSym=countLineFlag(file, "symmetries");
  printf("nsym=%d\n", nSym);
  
  Symmetries sym; allocateSymmetries(&sym, nSym, nSites); //kind of a constructor.
  //int symmetries[nSym*nSites];
  //Permutation *symmetries;
  //symmetries = (Permutation *) malloc(nSym * sizeof (Permutation *));
  readSymmetries(file, nSites, &sym);
  //sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
  printSymmetries(&sym); 
  
  return 0;
}


