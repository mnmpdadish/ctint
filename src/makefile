
MKL          = /opt/intel/mkl/lib/intel64
LIBS         = ${MKL}/libmkl_intel_lp64.a -Wl,--start-group $(MKL)/libmkl_blas95_lp64.a $(MKL)/libmkl_lapack95_lp64.a $(MKL)/libmkl_sequential.a ${MKL}/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 
#LIBS         = -llapack -lblas

all: build

build: 
	gcc -ffast-math  -Wall -O2 monte_carlo.c -lm  $(LIBS) -lgfortran -pg
	
	
	
	
