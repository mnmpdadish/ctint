MKL     = /opt/intel/mkl/lib/intel64
LIBS    = ${MKL}/libmkl_intel_lp64.a -Wl,--start-group $(MKL)/libmkl_blas95_lp64.a $(MKL)/libmkl_lapack95_lp64.a $(MKL)/libmkl_sequential.a ${MKL}/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 
#LIBS    = -llapack -lblas
INCLUDE =
OPTIONS = -Wall#-ffast-math#-pg

arrays=src/util/arrays

all: build 
build: 
	gcc $(OPTIONS) -o mc src/monte_carlo.c -lm $(LIBS) -lgfortran 

test:
	gcc $(OPTIONS) -o tests src/test.c -lm $(LIBS) -lgfortran

clean:
	rm mc tests

	
	
