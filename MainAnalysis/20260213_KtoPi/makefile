ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

ExecuteKtoPiAnalysis: KtoPiAnalysis.cpp
	 g++ -O3 -I. -I/data/yjlee/strangeness/Strangeness/CommonCode/include -I./include \
	    $(ROOTCFLAGS) \
	    KtoPiAnalysis.cpp \
	    /data/yjlee/strangeness/Strangeness/CommonCode/library/StrangenessMessenger.o \
	    -o ExecuteKtoPiAnalysis \
	    $(ROOTLIBS)
