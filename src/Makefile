# Compilation of MC-Integrated
CC = g++
OPTFLAGS = -O3 -fopenmp -mtune=native -flto -ffast-math -finline-functions -fprefetch-loop-arrays
CFLAGS = -g -Wall -Wextra -Wshadow -Wconversion -Wunreachable-code -mcmodel=medium -fmax-errors=1 -std=c++17
OBJECTS = MC.o main.o potentials.o energy.o moves.o print.o read.o sample.o
SOURCEFILES = MC.cpp potentials.cpp energy.cpp moves.cpp print.cpp read.cpp sample.cpp
MAINFILE = main.cpp
PROG = simulate

all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) -o $(PROG) $(CFLAGS) $(OPTFLAGS) $(OBJECTS)
MC.o: MC.cpp
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(SOURCEFILES)
main.o: main.cpp
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $(MAINFILE)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $<

install:
	mv $(PROG) ../
clean:
	rm -rf *.o $(PROG)

