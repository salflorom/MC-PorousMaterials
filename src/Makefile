# Compilation of MC-Integrated
CC = g++
CFLAGS = -g -Wall -Wextra -Wshadow -Wconversion -Wunreachable-code -O3 -mcmodel=medium -fmax-errors=1 -std=c++17
OBJECTS = MC.o main.o potentials.o energy.o moves.o print.o read.o sample.o
PROG = mcexecutable
all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) -o $(PROG) $(CFLAGS) $(OBJECTS)
MC.o: MC.cpp
	$(CC) -c $(CFLAGS) MC.cpp potentials.cpp moves.cpp print.cpp read.cpp sample.cpp
main.o: main.cpp
	$(CC) -c $(CFLAGS) main.cpp

%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

install:
	mv $(PROG) ../
clean:
	rm -rf *.o $(PROG)

