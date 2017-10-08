OPTIONS = -std=c99 -pedantic -Wall 
EXEC = ctINT

all: executable
executable: src/test.c
	gcc $(OPTIONS) -o $(EXEC) src/test.c -lm -llapack -lblas -lgfortran

