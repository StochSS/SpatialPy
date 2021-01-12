CC = g++
CFLAGS = -c -std=c++14 -Wall -O3 -g
CFILES = $(wildcard *.c)
CPPFILES = $(wildcard *.cpp)
OBJFILES = $(CFILES:.c=.o) $(CPPFILES:.c=.o)
SIMFLAGS = -l. -std=c++14 -Wall -O3
LFLAGS = -pthread -lm
DSFMTFLAGS = -DHAVE_SSE2 -DDSFMT_MEXP=521
INCDIR = $(ROOTINC)/include/ $(ROOTINC)/external/ $(ROOTINC)/external/ANN/include/
INCDIRPARAMS = $(INCDIR:%=-I%)
OBJ = particle.o simulate.o count_cores.o output.o simulate_rdme.o simulate_threads.o model.o pthread_barrier.o
#OBJ = linked_list.o particle.o simulate.o count_cores.o output.o simulate_rdme.o binheap.o simulate_threads.o model.o pthread_barrier.o


.PHONY: all

all: ssa_sdpd

main.o:
	$(CC) -c -o main.o $(MODEL) $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

dSFMT.o:
	$(CC) -c -o dSFMT.o $(ROOTINC)/external/dSFMT/dSFMT.c $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

ANN.o:
	$(CC) -c -o ANN.o $(ROOTINC)/external/ANN/src/ANN.cpp $(INCDIRPARAMS) $(CFLAGS)

%.o: $(ROOT)/src/%.cpp
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

%.o: $(ROOT)/src/%.c
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS) $(DSFMTFLAGS)

ssa_sdpd: ANN.o dSFMT.o main.o $(OBJ)
	$(CC) -o ssa_sdpd $(OBJ)  main.o dSFMT.o ANN.o $(LFLAGS)
