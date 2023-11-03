# Compilation of MC-Integrated
SRC=src/

all:
	cd $(SRC); $(MAKE)
install:
	cd $(SRC); $(MAKE) install; rm MC.o main.o
clean:
	cd $(SRC); $(MAKE) clean

