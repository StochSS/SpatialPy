CC = gcc
CFLAGS = -O3 -Wall
LFLAGS = -pthread -lm
DSFMTFLAGS = -DHAVE_SSE2 -DDSFMT_MEXP=521
INCDIR = $(ROOTINC)/include/ $(ROOTINC)/external/
INCDIRPARAMS = $(INCDIR:%=-I%)
OBJ = linked_list.o particle.o simulate.o count_cores.o output.o simulate_rdme.o binheap.o simulate_threads.o model.o pthread_barrier.o


.PHONY: all

all: ssa_sdpd

main.o:
	$(CC) -c -o main.o $(MODEL) $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

%.o: $(ROOT)/src/%.c
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

ssa_sdpd: main.o $(OBJ)
	$(CC)  -o ssa_sdpd $(OBJ) $(ROOTINC)/external/dSFMT/dSFMT.c main.o $(DSFMTFLAGS) $(LFLAGS)
