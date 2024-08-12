
SOURCES = matmul.f90 time.f90
FILE ?= main.f90
TARGET ?= matmul
FFLAGS = -O3
FC ?= gfortran

default: run

clean:
	rm -f *.mod $(TARGET)

build: clean
	$(FC) $(FFLAGS) $(SOURCES) $(FILE) -o $(TARGET)

run: build
	./$(TARGET) 1500
