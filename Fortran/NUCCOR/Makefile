FC = gfortran -O3 -fopenmp
objects = mtc_patch.o mtc.o mtc_openmp.o mtc_main.o

INCLUDE=

all: ${objects}
	${FC} -o mtc.x ${objects}

%.o: %.f90
	${FC} -c $< ${INCLUDE}

%.o: %.F90
	${FC} -c $< ${INCLUDE}

.PHONY. : clean
clean:
	rm *.o *.mod mtc.x

