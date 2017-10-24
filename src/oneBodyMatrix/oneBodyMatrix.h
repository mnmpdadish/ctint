#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../util/stringUtil.h"
#include "../util/arrays/array.h"
#include "../matrixDouble.h"
#include "../matrixComplex.h"


typedef struct {
  int x;
  int y;
  int z;
} IntPosition;

// c = a+f*b
IntPosition addIntPosition(IntPosition a, IntPosition b, int f) {
  IntPosition c;
  c.x = a.x + f*b.x;
  c.y = a.y + f*b.y;
  c.z = a.z + f*b.z;
  return c;
}

int isEqual_IntPosition(IntPosition a, IntPosition b) {
  if(a.x != b.x) return 0;
  if(a.y != b.y) return 0;
  if(a.z != b.z) return 0;
  return 1;
}

int isZero_IntPosition(IntPosition a) {
  if(a.x != 0) return 0;
  if(a.y != 0) return 0;
  if(a.z != 0) return 0;
  return 1;
}

typedef struct {
  unsigned int n;
  unsigned int capacity;
  IntPosition * positions;
} MultiplePositions;


int init_MultiplePositions(MultiplePositions * mPositions) {
  mPositions->capacity = INIT_CAPACITY;
  mPositions->positions = malloc(mPositions->capacity * sizeof (IntPosition));
  mPositions->n=0;
  return 0;
}

void addSites_MultiplePositions(MultiplePositions * mPositions, int pos3D[3]) {
  unsigned int n=mPositions->n;
  if(mPositions->n > mPositions->capacity) 
    mPositions->positions = realloc(mPositions->positions, (mPositions->capacity *= 2) * sizeof(IntPosition));
  mPositions->positions[n].x=pos3D[0];
  mPositions->positions[n].y=pos3D[1];
  mPositions->positions[n].z=pos3D[2];
  mPositions->n++;  
}

  

void readSites(FILE * file, MultiplePositions * mPositions, char * flag){
  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0;
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        if(strBeginWithToken(tempbuff,flag)) found=1; 
      }
      else{
        unsigned int nElement = countElementInStr(tempbuff, " \t\n");
        if(tempbuff[0] == '#') continue;
        else if(((tempbuff[0] != '*') || (tempbuff[0] != '\n')) && nElement != 0) {
          
          int pos3D[3];
          int nRead;
          nRead = sscanf(tempbuff,"(%d,%d,%d) \n", &pos3D[0],&pos3D[1],&pos3D[2]);
          addSites_MultiplePositions(mPositions, pos3D);
          
          if(nRead!=3) {
            printf("Cannot read correctly the %s line: \n%s", flag, tempbuff); 
            exit(1);
          }
        }
        else break;
      }
    }
  }
}



void print_MultiplePositions(MultiplePositions * mPositions, char *flag) {
  unsigned int i;
  for(i=0;i<mPositions->n;i++) {
    printf("%s %d with position: (%d,%d,%d)\n", flag, i, mPositions->positions[i].x
                                                       , mPositions->positions[i].y
                                                       , mPositions->positions[i].z);
  }
}


void free_MultiplePositions(MultiplePositions * mPositions) {
  free(mPositions->positions);
}








typedef struct {
  IntPosition R; //cluster origine position.
  IntPosition r; //site position (relative to cluster origine)
} Folding;


int tripleProductInteger(IntPosition a, IntPosition b, IntPosition c) {
  return a.x*(b.y*c.z - c.y*b.z) - a.y*(b.x*c.z - c.x*b.z) + a.z*(b.x*c.y - b.y*c.x);
} 

Folding Fold(IntPosition position, MultiplePositions sites, MultiplePositions superlattice){
  assert(superlattice.n==3);
  int volume = tripleProductInteger(superlattice.positions[0], superlattice.positions[1], superlattice.positions[2]);
  if(volume<0) {
    volume *=-1;
    superlattice.positions[0].x*=-1; 
    superlattice.positions[0].y*=-1; 
    superlattice.positions[0].z*=-1;
  }
  int found = 0, i ;
  IntPosition R, Q;//, Zeros;  
  //Zeros.x=0; Zeros.y=0; Zeros.z=0;
  for(i=0; i<sites.n; i++){
    IntPosition r = addIntPosition(position, sites.positions[i] ,-1);
    R.x = tripleProductInteger(r, superlattice.positions[1], superlattice.positions[2]);
    R.y = tripleProductInteger(superlattice.positions[0], r, superlattice.positions[2]);
    R.z = tripleProductInteger(superlattice.positions[0], superlattice.positions[1], r);
    Q.x = R.x % volume;
    Q.y = R.y % volume;
    Q.z = R.z % volume;
    if(isZero_IntPosition(Q) ){
      found=1;
      break;
    }
  }
  if(found){
    IntPosition foldingR;
    foldingR.x = (superlattice.positions[0].x*R.x + superlattice.positions[1].x*R.y + superlattice.positions[2].x*R.z)/volume;
    foldingR.y = (superlattice.positions[0].y*R.x + superlattice.positions[1].y*R.y + superlattice.positions[2].y*R.z)/volume;
    foldingR.z = (superlattice.positions[0].z*R.x + superlattice.positions[1].z*R.y + superlattice.positions[2].z*R.z)/volume;
    
    Folding foldingResult;
    foldingResult.R = foldingR; 
    foldingResult.r=sites.positions[i];
    
    return foldingResult;
  }
  else{
    printf("error: Unable to fold the position (%d,%d,%d) on the original cluster with this superlattice", position.x, position.y, position.z);
    exit(1);
  }
  
  
}


/*  
   vol = (dot(e[0],cross(e[1],e[2])))   # equivalent to determinant, but keep the number integers
   if vol < 0:
      vol *= -1
      e[0] *= -1 
   e = transpose(e)
   found = False
   
   for siteIndex in range(0,len(sites)):
      site = sites[siteIndex]
      r = position - site
      R = CramersRuleInt(e,r)
      Q = R % vol
      if all(Q == array([0,0,0])):
         found= True
         break   
   if found:
      R2 = dot(e,R)//vol
      return R2,site
   else:
      print 'Unable to fold the position',position, 'on the original cluster with this superlattice'
      assert()

*/

// ------------------------------------------------------------------

typedef struct {
  unsigned int *index1;
  unsigned int *index2;
  double *value;
  IntPosition *clusterPosition;
  unsigned int n;
  unsigned int capacity;
} Sparse_tMatrix; // contains a list of 

// ------------------------------------------------------------------

typedef struct {
  char label[64]; //name should be under 64 characters.
  double coefficient;
  IntPosition posDiff; 
} OperatorDef;


typedef struct {
  unsigned int n;
  unsigned int capacity;
  unsigned int nSites;
  OperatorDef * operators;
  Sparse_tMatrix sparse;
} tMatrix;


int init_tMatrix(tMatrix * tMat) {
  tMat->capacity = INIT_CAPACITY;
  tMat->operators = (OperatorDef *) malloc(tMat->capacity * sizeof (OperatorDef ));
  tMat->n=0;
  
  tMat->sparse.capacity = INIT_CAPACITY;
  tMat->sparse.n = 0;
  tMat->sparse.index1 = (unsigned int *) malloc(tMat->sparse.capacity * sizeof (unsigned int));
  tMat->sparse.index2 = (unsigned int *) malloc(tMat->sparse.capacity * sizeof (unsigned int));
  tMat->sparse.value  = (double *)       malloc(tMat->sparse.capacity * sizeof (double));
  tMat->sparse.clusterPosition = (IntPosition *) malloc(tMat->sparse.capacity * sizeof (IntPosition));
  return 0;
}

void addOperators_tMatrix(tMatrix * tMat, char label[64], double coefficient, int posDiff3D[3]) {
  if(tMat->n > tMat->capacity) tMat->operators = realloc(tMat->operators, (tMat->capacity *= 2) * sizeof(OperatorDef));
  memset(tMat->operators[tMat->n].label,0,64);
  strcpy(tMat->operators[tMat->n].label,label);
  tMat->operators[tMat->n].coefficient=coefficient;
  tMat->operators[tMat->n].posDiff.x=posDiff3D[0];
  tMat->operators[tMat->n].posDiff.y=posDiff3D[1];
  tMat->operators[tMat->n].posDiff.z=posDiff3D[2];
  tMat->n++;  
}

  

void readOperators_tMatrix(FILE * file, tMatrix * tMat){
  rewind(file);
  char tempbuff[256];  //each line should not be above 256 char long.
  int found=0;
        
  while(!feof(file)) 
  {
    if (fgets(tempbuff,256,file)) {
      if(found==0){
        if(strBeginWithToken(tempbuff,"one-body")) found=1; 
      }
      else{
        unsigned int nElement = countElementInStr(tempbuff, " \t\n");
        if(tempbuff[0] == '#') continue;
        else if(((tempbuff[0] != '*') || (tempbuff[0] != '\n')) && nElement != 0) {
          
          char label[64];
          float coefficient;
          int posDiff3D[3];
          int minusPosDiff3D[3];
          int nRead,i;
          nRead = sscanf(tempbuff,"%s (%d,%d,%d) %f \n",label, &posDiff3D[0],&posDiff3D[1],&posDiff3D[2],&coefficient);
          addOperators_tMatrix(tMat, label, coefficient, posDiff3D);
          int isNot_0_0_0=0;
          for(i=0;i<3;i++) {
            minusPosDiff3D[i]=-posDiff3D[i];
            if(minusPosDiff3D[i] != posDiff3D[i]) isNot_0_0_0=1;
          }
          if(isNot_0_0_0) addOperators_tMatrix(tMat, label, coefficient, minusPosDiff3D);
          
          if(nRead!=5) {
            printf("Cannot read correctly the one-body line: \n%s", tempbuff); 
            exit(1);
          }
        }
        else break;
      }
    }
  }
}

void addElementSparse_tMatrix(tMatrix * tMat, unsigned int i1, unsigned int i2, IntPosition foldingR, double value) {
  unsigned int N = tMat->sparse.n;
  if(N == tMat->sparse.capacity) {
    tMat->sparse.capacity *= 2;
    //printf("new capacity = %d",tMat->sparse.capacity);
    tMat->sparse.index1 = realloc(tMat->sparse.index1, (tMat->sparse.capacity ) * sizeof(unsigned int));
    tMat->sparse.index2 = realloc(tMat->sparse.index2, (tMat->sparse.capacity) * sizeof(unsigned int));
    tMat->sparse.value  = realloc(tMat->sparse.value, (tMat->sparse.capacity) * sizeof(double));
    tMat->sparse.clusterPosition = realloc(tMat->sparse.clusterPosition, (tMat->sparse.capacity) * sizeof(IntPosition));
  }
  tMat->sparse.index1[N] = i1;
  tMat->sparse.index2[N] = i2;
  tMat->sparse.clusterPosition[N] = foldingR;
  tMat->sparse.value[N] = value;
  tMat->sparse.n++;
}

void defineSparse_tMatrix(tMatrix * tMat, MultiplePositions *sites, MultiplePositions *superlattice) {
  unsigned int i,j,k;
  tMat->nSites = sites->n;
  for(i=0;i<tMat->n;i++){
    for(j=0;j<sites->n;j++){
      IntPosition newPosition1 = addIntPosition(sites->positions[j], tMat->operators[i].posDiff, 1);
      Folding folding = Fold(newPosition1, *sites, *superlattice);
      unsigned int newIndex=-1;
      for(k=0;k<sites->n;k++) if(isZero_IntPosition(addIntPosition(sites->positions[k],  folding.r, -1))) newIndex = k;
      addElementSparse_tMatrix(tMat,newIndex,j,folding.R,tMat->operators[i].coefficient); 
    }
  }
}

void print_tMatrix(tMatrix * tMat) {
  unsigned int i;
  for(i=0;i<tMat->n;i++) {
    printf("operator %d: '%s' with coeff %f and position difference vector: (%d,%d,%d)\n", 
            i, tMat->operators[i].label, tMat->operators[i].coefficient,
            tMat->operators[i].posDiff.x,
            tMat->operators[i].posDiff.y,
            tMat->operators[i].posDiff.z);
  }
  printf("\n");
  for(i=0;i<tMat->sparse.n;i++) {
    printf("sparse(%d,%d) = % 6.3f and cluster position: (%d,%d,%d)\n", 
            tMat->sparse.index1[i], tMat->sparse.index2[i], tMat->sparse.value[i],
            tMat->sparse.clusterPosition[i].x,
            tMat->sparse.clusterPosition[i].y,
            tMat->sparse.clusterPosition[i].z);
  }
  printf("\n");
}



void calculate_tMatrixK_2D(tMatrix * tMat, cMatrix * tMatrixK, double kx, double ky) {
  unsigned int i;
  assert(tMatrixK->N==tMat->nSites);
  reset_cMatrix(tMatrixK);
  for(i=0;i<tMat->sparse.n;i++){
    double phase = kx*tMat->sparse.clusterPosition[i].x+ky*tMat->sparse.clusterPosition[i].y;
    double complex ek = cos(phase)+I*sin(phase); 
  
    ELEM(tMatrixK, tMat->sparse.index1[i], tMat->sparse.index2[i]) += tMat->sparse.value[i] * ek;
  }
}


void calculate_hybFirstMoments(tMatrix * tMat, dMatrix * hybFM) {
  unsigned int i,j,k,n1,n2; // i,j: the indices of the matrix, k is the dummy summation index
  assert(hybFM->N==tMat->nSites);
  reset_dMatrix(hybFM);
  for(i=0;i<tMat->sparse.n;i++)
    for(j=0;j<tMat->sparse.n;j++)
      for(k=0;k<tMat->sparse.n;k++)
        for(n1=0;n1<tMat->sparse.n;n1++) 
          if( tMat->sparse.index1[n1]==i && tMat->sparse.index2[n1]==k && !isZero_IntPosition(tMat->sparse.clusterPosition[n1]) ) {
            for(n2=0;n2<tMat->sparse.n;n2++)
              if( tMat->sparse.index1[n2]==k && tMat->sparse.index2[n2]==j && !isZero_IntPosition(tMat->sparse.clusterPosition[n2]) )
                if(isZero_IntPosition(addIntPosition(tMat->sparse.clusterPosition[n1],  tMat->sparse.clusterPosition[n2], 1)))
                  ELEM(hybFM, i,j) += tMat->sparse.value[n1] * tMat->sparse.value[n2];
          }
}


/*
      for i1 in range(0,self.n_label):
         for i2 in range(0,self.n_label):
            for dum in range(0,self.n_label):
               keys1 = (i1,dum) 
               keys2 = (dum,i2)
               if (keys1 in sparseK) and (keys2 in sparseK):
                  for nn1 in range(0,len(sparseK[keys1])):
                     for nn2 in range(0,len(sparseK[keys2])):
                        if all((sparseK[keys1][nn1][1]+sparseK[keys2][nn2][1])==array([0,0,0])) :
                           matrixHybFirstMoment[i1,i2] += ((sparseK[keys1])[nn1][0])*((sparseK[keys2])[nn2][0])
      


      #precalculate the k-indep matrix_0
      matrix_0 =  zeros((self.n_label,self.n_label),dtype=complex)
      for keys in self.sparse_0:
         matrix_0[keys[0],keys[1]]=self.sparse_0[keys]
         matrix_0[keys[1],keys[0]]=self.sparse_0[keys]
      self.matrix_0 = matrix_0
*/


void free_tMatrix(tMatrix * tMat) {
  free(tMat->operators);
  free(tMat->sparse.index1);
  free(tMat->sparse.index2);
  free(tMat->sparse.value);
  free(tMat->sparse.clusterPosition);
}


