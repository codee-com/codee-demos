FC ?= gfortran
FFLAGS ?= -O3
SRC = himeno.f90
TARGET ?= himeno

all: clean build

build: $(SRC)
	$(FC) -o $(TARGET) $^ $(FFLAGS)

run: build
	./$(TARGET)

clean:
	rm -f bound.mod  bound.pw.mod  mtrx.mod  mtrx.pw.mod  others.mod  others.pw.mod  pres.mod  pres.pw.mod  work.mod  work.pw.mod

