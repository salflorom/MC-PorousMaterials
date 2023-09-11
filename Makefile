# Compilation of MC-Integrated
SRC=src/

all:
	cd $(SRC); $(MAKE)
install:
	cd $(SRC); $(MAKE) install
clean:
	cd $(SRC); $(MAKE) clean

