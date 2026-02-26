

default: ExecuteClosure
	./ExecuteClosure --Fraction 0.01

Execute: DeriveEfficiency.cpp
	g++ DeriveEfficiency.cpp -o Execute `root-config --cflags --libs` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/StrangenessMessenger.o

ExecuteClosure: ClosureCheck.cpp
	g++ ClosureCheck.cpp -o ExecuteClosure `root-config --cflags --libs` \
		-I$(ProjectBase)/CommonCode/include $(ProjectBase)/CommonCode/library/StrangenessMessenger.o
