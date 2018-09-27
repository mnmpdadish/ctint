#!/bin/bash

mc_DIR=../../
MODEL=dmft

ITER=0
ITERMAX=100

if [ -a logfile ]
  then rm logfile
fi

while [ $ITER -le $ITERMAX ]
do
  echo begin iteration $ITER at: `date` >> logfile 
  ${mc_DIR}/mc ${MODEL}.model ${ITER}
  
  ITER=$[$ITER+1]
done

