CC=gcc -fopenmp
CFLAGS=-O3 -Wall
LDFLAGS=-lm

OBJS=qopt.o
BIN=qopt

$(BIN): $(OBJS)
	$(CC) $^ $(LDFLAGS) -o $(BIN) 

mpi: MPI-qopt.c qopt-lib.o
	mpicc $(CFLAGS) $^ $(LDFLAGS) -o MPI-qopt

clean: 
	rm $(BIN) $(OBJS)
