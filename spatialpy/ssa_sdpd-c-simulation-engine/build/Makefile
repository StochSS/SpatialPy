CC = g++
CFLAGS = -c -std=c++14 -Wall -O3 -g
SIMFLAGS = -l. -std=c++14 -Wall -O3
LFLAGS = -pthread -lm
INCDIR = $(ROOTINC)/include/ $(ROOTINC)/external/ $(ROOTINC)/external/ANN/include/
INCDIRPARAMS = $(INCDIR:%=-I%)
OBJ = particle.o simulate.o count_cores.o output.o simulate_rdme.o simulate_threads.o model.o pthread_barrier.o
#OBJ = linked_list.o particle.o simulate.o count_cores.o output.o simulate_rdme.o binheap.o simulate_threads.o model.o pthread_barrier.o


.PHONY: all

all: ssa_sdpd

main.o:
	$(CC) -c -o main.o $(MODEL) $(INCDIRPARAMS) $(CFLAGS)

ANN.o:
	$(CC) -c -o ANN.o $(ROOTINC)/external/ANN/src/ANN.cpp $(INCDIRPARAMS) $(CFLAGS)

%.o: $(ROOT)/src/%.cpp
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS)

%.o: $(ROOT)/src/%.c
	$(CC) -c -o $@ "$<" $(INCDIRPARAMS) $(CFLAGS)

ssa_sdpd: ANN.o main.o $(OBJ)
	$(CC) -o ssa_sdpd $(OBJ)  main.o ANN.o $(LFLAGS)
