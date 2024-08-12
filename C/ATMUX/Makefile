SOURCES = lib/Matrix2D.c lib/Vector.c lib/CRSMatrix.c
FILE ?= atmux.c
TARGET ?= atmux
CFLAGS = -std=c99 -O3 -Ilib -fopenmp

default: run

clean:
	rm -f $(TARGET)

build: clean
	$(CC) $(FILE) $(SOURCES) $(CFLAGS) -o $(TARGET)

run: build
	./$(TARGET) 10000
	