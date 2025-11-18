FC = gfortran
CC = gcc

FFLAGS ?= -c -cpp -Wall -Wuninitialized -Wmaybe-uninitialized
CFLAGS ?= -c -Wall -Wuninitialized -Wmaybe-uninitialized

F90_SRCS := $(wildcard *.f90)
C_SRCS  := $(wildcard *.c)

PROGS := $(F90_SRCS:.f90=) $(C_SRCS:.c=)

.PHONY: all clean

all: $(PROGS)

%: %.f90
	$(FC) $(FFLAGS) $< -o $@

%: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm -f $(PROGS)
