#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../util/stringUtil.h"
#include "../util/arrays/array.h"


typedef struct {
  char label;
  double coefficient;
  int positionVector[3]; 
} Operator;


typedef struct {
  unsigned int n;
  unsigned int capacity;
  Operator * operators;
  //double complex * 
} tMatrix;


