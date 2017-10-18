#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "../util/stringUtil.h"
#include "../util/arrays/array.h"


typedef struct {
  char name;
  double coefficient;
  int positionVector[3]; 
} Operator;


typedef struct {
  unsigned int n;
  unsigned int capacity;
  Operator * operators;
} tMatrix;


