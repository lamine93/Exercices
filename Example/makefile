PROGRAM = vecmat_gpun.x
CC = gcc
FC = gfortran
LIBS = 
LDR = ${FC} 
FFLAGS = -fopenmp -O3 -c 
CFLAGS=  -fopenmp 
LDFLAGS= -fopenmp


OBJS = vecmat_gpun.o

.f.o:
	$(FC) $(FFLAGS) $*.f
.c.o:
	$(CC) $(CFLAGS) $*.c

all:	$(PROGRAM)

$(PROGRAM):	$(OBJS) 
	$(LDR) $(LDFLAGS) -o $(PROGRAM)  $(OBJS) $(LIBS)



