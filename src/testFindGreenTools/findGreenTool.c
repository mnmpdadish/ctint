#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "stringUtil.h"



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
  assert(nSites<=64);
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
    assert(perm->data[i] != value );
  }
  //printf("%d %d\n", value, perm->nSites);
  assert(value<perm->nSites && value>=0);

  perm->data[perm->nAssigned]=value;
  perm->nAssigned++;
  assert(perm->nAssigned <= perm->nSites);
  return 0;
}


char charHexa(int digit) {
  //any cluster should not have more than 64 sites, even in 2050 xD
  char characters[] = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ!?"; 
  char char1;
  if(digit < 64 && digit >= 0) char1 = characters[digit];
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






/*
int countNumberOfElementInOneLine(const char *str1){
  char str2[256]="";
  strncpy(str2, str1, 256); //do a copy
  char *token, *saveptr;
  const char delimiter = ' ';
  token = strtok_r(str2, &delimiter, &saveptr);

  int N=0;
  while( token != NULL && *token != '\n') {
    N++;
    token = strtok_r(NULL, &delimiter, &saveptr);
  }
  return N;
}
*/


int countLineFlag(FILE * file, char *flag) {

  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, N=0;
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        char tmpstr1[50]; tmpstr1[0]='_';
        sscanf(tempbuff, "%256s\n", tmpstr1);
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
    //printf("\n");
  }
  return N;
}


void readSymmetries(FILE * file, int nSites, Symmetries *sym) {

  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0, i=0;
  char name[]="symmetries";
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      //printf("%d seen= %s",found, tempbuff);
      if(found==0){
        char tmpstr1[50]; tmpstr1[0]='_';
        //printf("word= %s   fisrtletter=%c\n", tmpstr1, tempbuff[0]);
        //exit(1);
        sscanf(tempbuff, "%s\n", tmpstr1);
        //printf("word= %s   fisrtletter=\n", tmpstr1, tempbuff[0]);
        if (strcmp(tmpstr1,name)==0) { 
          found=1;
          continue;
        }
      }
      else{
        if(tempbuff[0] == '#') continue;
        else if((tempbuff[0] != '\n') && countNumberOfElementInOneLine(tempbuff) != 0) {
		      unsigned int nElement = countNumberOfElementInOneLine(tempbuff);
          int arrayInt[nElement];
          readIntInOneLine(tempbuff, arrayInt);
          int j;
          for(j=0;j<nElement;j++) {
            //printf("number= %d\n",arrayInt[j]);
            addElementToPermutation(&(sym->permutation[i]),arrayInt[j]);
          }
				    
          
          
		      //else break;
		      /*
		      printf("tempbuff= %s",tempbuff);
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
				  */
				  //printf("tempbuff= %s",tempbuff);
				  //printf("\n %d %d %d\n",sym->permutation[i].nSites,sym->n,nElement);
				  assert(sym->permutation[i].nSites == nElement);
				  i++;
		    }
		    else break;
        //else if((tempbuff[0] != '\n') && countNumberOfElementInOneLine(tempbuff) != 0) N++;
        //else break;
      }
    }
  }   
}


typedef struct {
  unsigned int *i;
  unsigned int *j;
  unsigned int nElement;
  unsigned int nSites;
} GreenSym;

int buildGreenSym(GreenSym * greenSym, int nSites) {
  greenSym->i = (unsigned int *) malloc(nSites*nSites * sizeof (unsigned int));
  greenSym->j = (unsigned int *) malloc(nSites*nSites * sizeof (unsigned int));
  greenSym->nElement = nSites*nSites;
  greenSym->nSites = nSites;
  
  int i,j;
  for(i=0;i<nSites;i++){
    for(j=0;j<nSites;j++){
      greenSym->i[i*nSites+j]=(i<j)? i:j;
      greenSym->j[i*nSites+j]=(i<j)? j:i;
    }
  } 
  return 0;
}

int printGreenSym(GreenSym * greenSym) {  
  int i,j,N=greenSym->nSites;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      printf("%c%c ",charHexa(greenSym->i[N*i+j]),charHexa(greenSym->j[N*i+j]));
    }
    printf("\n");
  } 
  return 0;
}


int countIndependantGreen(GreenSym * greenSym) {  
  int i,j,N=greenSym->nSites, nIndep=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(greenSym->i[N*i+j]==i && greenSym->j[N*i+j]==j) nIndep++;
    }
  } 
  return nIndep;
}



int symmetrizeOneGreenELement(GreenSym * greenSym, Symmetries * sym) {  
  int i,j,k, N=greenSym->nSites;
  for(k=0; k<sym->n; k++){
    for(i=0;i<greenSym->nSites;i++){
      for(j=0;j<greenSym->nSites;j++){
        int permutation_i = i;
        int permutation_j = j;
        int isFirst=1, n=0;
        while( ((permutation_i != i) || (permutation_j != j)) || isFirst) { //while didn't loop on orbit
          n++;
          permutation_i= sym->permutation[k].data[permutation_i];
          permutation_j= sym->permutation[k].data[permutation_j];
          int refIndex=N*i+j ;
          int newIndex=N*permutation_i+permutation_j;
          if( (greenSym->i[newIndex] != greenSym->i[refIndex]) || (greenSym->j[newIndex] != greenSym->j[refIndex]) ) {
            //printf("k=%d  i=%d j=%d %d %d  %c%c %c%c  %d %d\n",k,i,j,permutation_i,permutation_j, charHexa(greenSym->i[refIndex]), charHexa(greenSym->j[refIndex]), charHexa(greenSym->i[newIndex]), charHexa(greenSym->j[newIndex]), refIndex, newIndex);
            if( N*greenSym->i[refIndex] + greenSym->j[refIndex] < N*greenSym->i[newIndex] + greenSym->j[newIndex] ) {
            //refIndex < newIndex){
              greenSym->i[newIndex]=greenSym->i[refIndex];
              greenSym->j[newIndex]=greenSym->j[refIndex];
            }
            else{
              greenSym->i[refIndex]=greenSym->i[newIndex];
              greenSym->j[refIndex]=greenSym->j[newIndex];
            }
            //printf("k=%d  i=%d j=%d %d %d  %c%c %c%c\n",k,i,j,permutation_i,permutation_j, charHexa(greenSym->i[refIndex]), charHexa(greenSym->j[refIndex]), charHexa(greenSym->i[newIndex]), charHexa(greenSym->j[newIndex]));

            return 1;
          }
          isFirst=0;
          if(n>20) exit(1);
        }
      }
    }
  }
  return 0;
}


int main() {
  int nSites=64;
  char fileName[]="plaquette8x8.in";
  FILE * file = fopen(fileName, "rt");
  if(file == NULL) {printf("file %s not found", fileName); exit(1);}
  printf("reading symmetries from %s:\n", fileName);

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
  
  GreenSym greenSym;
  buildGreenSym(&greenSym,nSites);
  printf("\n");
  printGreenSym(&greenSym);
  while(symmetrizeOneGreenELement(&greenSym,&sym)){}
  printf("\n");
  printGreenSym(&greenSym);
  
  printf("nIndep=%d\n",countIndependantGreen(&greenSym));

  return 0;
}


