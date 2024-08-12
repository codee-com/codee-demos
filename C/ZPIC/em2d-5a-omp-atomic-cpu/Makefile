## GCC options
CC = gcc
CFLAGS = -O3 -std=c99 -pedantic -fopenmp

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic /Qopenmp

# Clang options
#CC = clang
#CFLAGS = -Ofast -fopenmp

# Portland options
#CC = pgcc
#CFLAGS = -Ofast -mp

# Cray options
#CC = craycc
#CFLAGS = -Ofast -h omp, noacc

LDFLAGS = -lm


SOURCE = current.c emf.c particles.c random.c timer.c main.c simulation.c zdf.c one-malloc.c

TARGET = zpic

OBJ = $(SOURCE:.c=.o)

all : $(SOURCE) $(TARGET)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $@

.c.o:
	$(CC) -c $(CFLAGS) $< $(LDFLAGS) -o $@

clean:
	@touch $(TARGET) $(OBJ)
	rm -f $(TARGET) $(OBJ)
