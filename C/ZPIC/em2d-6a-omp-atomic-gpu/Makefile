# GCC options
#CC = gcc
#FLAGS = -Ofast -std=c99 -pedantic -fopenmp -foffload=nvptx-none="-Ofast -lm -misa=sm_35"

# Intel icc compiler
#CC = icc
#CFLAGS = -restrict -Ofast -std=c99 -pedantic /Qopenmp

# Clang options
CC = clang
CFLAGS = -Ofast -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda

# Portland options
#CC = pgcc
#CFLAGS = -Ofast -mp

# Cray options
#CC = craycc
#CFLAGS = -Ofast -h omp, noacc

# NVC options
#CC = nvc
#CFLAGS = -fast -Minfo=mp -mp -gpu=cc70,cuda11.0

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
