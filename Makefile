# Compilation of MC-Integrated
CC = g++
CFLAGS = -g -Wall -Wextra -Wshadow -Wconversion -Wunreachable-code -O3 -mcmodel=medium
OBJECTS = MC.o main.o
PROG = mcexecutable
all: $(PROG)

$(PROG): $(OBJECTS)
	$(CC) -o $(PROG) $(CFLAGS) $(OBJECTS)
MC.o: MC.cpp
	$(CC) -c $(CFLAGS) MC.cpp
main.o: main.cpp
	$(CC) -c $(CFLAGS) main.cpp

%.o: %.cpp
	$(CC) -c $(CFLAGS) $<

clean:
	rm -rf *.o $(PROG)

