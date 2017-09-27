OPTIONS = -Wall
EXEC = ctINT

all: executable
executable: src/main.c
	gcc $(OPTIONS) -o $(EXEC) src/matrix.c -llapack -lblas -lm -lgfortran -lstdc++ 

