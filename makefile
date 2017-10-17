EXEC = ctINT
INCLUDE =
OPTIONS = -O2 -std=c99 -pedantic -Wall $(INCLUDE)

arrays=src/util/arrays

all: build findGreenTool
build: 
	gcc $(OPTIONS) -o tests src/test.c -lm -llapack -lblas -lgfortran

findGreenTool: 
	gcc $(OPTIONS) -o findGreenTool src/findGreenTools/findGreenTool.c 

clean:
	rm findGreenTool tests
