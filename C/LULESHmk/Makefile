SHELL=/bin/sh
CLASS=W
VERSION=
SFILE=config/suite.def

CC = gcc
CFLAGS = -O2 -fopenmp

CPP =  src/luleshmk.c
OBJ = src/luleshmk

all: clean main

main: 
	$(CC) -O2 $(CPP) -o $(OBJ) -lm

run:
	./$(OBJ)

# Awk script courtesy cmg@cray.com, modified by Haoqiang Jin
suite:
	@ awk -f sys/suite.awk SMAKE=$(MAKE) $(SFILE) | $(SHELL)

clean:
	rm -f $(OBJ)