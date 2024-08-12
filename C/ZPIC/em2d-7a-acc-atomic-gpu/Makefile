# GCC options
#CC = gcc
#CFLAGS = -O3 -std=c99 -pedantic -fopenacc -foffload=nvptx-none="-lm"

# Cray options
#CC = craycc
#CFLAGS = -Ofast -h acc, noomp

# PGI options
#CC = pgcc
#CFLAGS = -Ofast -acc

# NVC options
CC = nvc
CFLAGS = -fast -Minfo=acc -acc -gpu=cc70,cuda11.0

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
