SHELL=/bin/sh
BENCHMARK=cg
BENCHMARKU=CG
FILE ?= cg.c
TARGET ?= cg

include ../config/make.def

include ../sys/make.common

OBJS = $(TARGET).o \
       ${COMMON}/print_results.o  \
       ${COMMON}/${RAND}.o \
       ${COMMON}/c_timers.o \
       ${COMMON}/wtime.o


${PROGRAM}: config ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${PROGRAM} ${OBJS} ${C_LIB}

.c.o:
	${CCOMPILE} $<

$(TARGET).o:		$(FILE)  globals.h npbparams.h

clean:
	- rm -f *.o *~ 
	- rm -f npbparams.h core
