
#LIBS    =-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
MKL     = /opt/intel/mkl/lib/intel64
LIBS    = ${MKL}/libmkl_intel_lp64.a -Wl,--start-group $(MKL)/libmkl_blas95_lp64.a $(MKL)/libmkl_lapack95_lp64.a $(MKL)/libmkl_sequential.a ${MKL}/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 
#LIBS    = -llapack -lblas
#LIBS    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl
#LIBS    = ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 

#LIBS    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -mkl -lpthread -ldl
INCLUDE =
OPTIONS = -Wall -O3 #-ffast-math#-pg

arrays=src/util/arrays

all: build 
build: 
	gcc $(OPTIONS) -o mc src/monte_carlo.c -lm $(LIBS) -lgfortran 

test:
	gcc $(OPTIONS) -o tests src/test.c -lm $(LIBS) -lgfortran

clean:
	rm mc tests

	
	
