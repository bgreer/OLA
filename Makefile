
CC=mpicc -ip -ipo -O3
#CC=mpicc -O2 -g

MKL = -L/central/intel/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
FITS=-L/shared/cfitsio-3.24/lib -lcfitsio

ola : main.o kernels.o measurements.o function.o optimization.o interp.o parse_input.o
	$(CC) main.o kernels.o measurements.o function.o optimization.o interp.o parse_input.o -o ola $(FITS) -lm $(MKL) -lrt

main.o : main.c header.h
	$(CC) -c main.c

kernels.o : kernels.c header.h
	$(CC) -c kernels.c

measurements.o : measurements.c header.h
	$(CC) -c measurements.c

function.o : function.c header.h
	$(CC) -c function.c

optimization.o : optimization.c header.h
	$(CC) -c optimization.c

interp.o : interp.c header.h
	$(CC) -c interp.c

parse_input.o : parse_input.c header.h
	$(CC) -c parse_input.c

clean : 
	rm *.o ola

matrixtest : matrixtest.c
	$(CC) -o matrixtest matrixtest.c $(MKL) -lm
