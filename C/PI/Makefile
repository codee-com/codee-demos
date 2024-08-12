FILE ?= pi.c
TARGET ?= pi
CFLAGS = -fopenmp -O3 -lm

default: run

clean:
	rm -f $(TARGET)

build: clean
	$(CC) $(FILE) $(CFLAGS) -o $(TARGET)

run: build
	./$(TARGET) 1000000000
