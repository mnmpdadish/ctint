
#LIBS    =-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl
#MKL     = /opt/intel/mkl/lib/intel64
#LIBS    = ${MKL}/libmkl_intel_lp64.a -Wl,--start-group $(MKL)/libmkl_blas95_lp64.a $(MKL)/libmkl_lapack95_lp64.a $(MKL)/libmkl_sequential.a ${MKL}/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 
#LIBS    = -llapack -lblas
#LIBS    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl
#LIBS    = ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lm -ldl 

#LIBS    = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -mkl -lpthread -ldl

cluster := $(shell echo $(CLUSTER))

ifeq ($(cluster),graham)
	LIBS     = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -mkl -lpthread -ldl
	COMPILER = icc
else ifdef MKLROOT 
	LIBS     = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl 
	COMPILER = gcc
else 
	LIBS     = -llapack -lblas -lrt
	COMPILER = gcc
endif

#LIBS     = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -mkl -lpthread -ldl
#COMPILER = icc
#LIBS     = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -ldl 
LIBS     = -llapack -lblas -lrt
COMPILER = gcc

OPTIONS = -Wall -O0 

all: build 
build: 
	@echo "Compiling for $(cluster)..."
	$(COMPILER) $(OPTIONS) -o mc src/monte_carlo.c -lm $(LIBS) -lgfortran 

test:
	$(COMPILER) $(OPTIONS) -o tests src/test.c -lm $(LIBS) -lgfortran

clean:
	rm mc 

	
	
