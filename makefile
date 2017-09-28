OPTIONS = -std=c99 -pedantic -Wall 
EXEC = ctINT

all: executable
executable: src/main.c
	gcc $(OPTIONS) -o $(EXEC) src/matrix.c -lm -llapack -lblas -lgfortran

