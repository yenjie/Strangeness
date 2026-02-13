efault: all

all: Setup library/StrangenessMessenger.o

Setup:
	mkdir -p library
	mkdir -p binary

library/StrangenessMessenger.o: source/StrangenessMessenger.cpp include/StrangenessMessenger.h
	g++ source/StrangenessMessenger.cpp -Iinclude -c -o library/StrangenessMessenger.o `root-config --cflags`