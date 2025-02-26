FC ?= gfortran
FLAGS ?= -O3 -fopenmp
objects = mtc_patch.o mtc.o mtc_openmp.o mtc_main.o

INCLUDE=

all: ${objects}
	${FC} ${FLAGS} -o mtc.x ${objects}

%.o: %.f90
	${FC} ${FLAGS} -c $< ${INCLUDE}

%.o: %.F90
	${FC} ${FLAGS} -c $< ${INCLUDE}

.PHONY. : clean
clean:
	rm *.o *.mod mtc.x

