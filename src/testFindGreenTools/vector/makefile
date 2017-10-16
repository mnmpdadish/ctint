CC=gcc
CFLAGS= -O2
RM=rm -rf
OUT=testArray

all: build

build: main.o array.o
	$(CC) $(CFLAGS) -o $(OUT) test.c array.c
	$(RM) *.o

debug: CFLAGS+=-DDEBUG_ON
debug: build

main.o: test.c array.h
	$(CC) $(CFLAGS) -c test.c

vector.o: array.c array.h
	$(CC) $(CFLAGS) -c array.c

clean:
	$(RM) *.o $(OUT)
