OPTIONS = -O2 -Wall
EXEC = ctINT

all: executable
executable: src/main.c
	gcc $(OPTIONS) -o $(EXEC) src/main.c -lm  

