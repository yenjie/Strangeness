ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --glibs)

default: ExecuteMakeKStarSignalOnlyHistograms ExecuteFitKStarSignalOnlyShapes ExecuteMakeKStarWrongAssignmentHistograms ExecuteCombineKStarAssignments ExecuteMakeKStarSBHistograms ExecuteFitKStarSB ExecuteFitKStarSBWithCrossFeed ExecuteMakePhiWrongAsKStarHistograms ExecuteMakeKShortWrongAsKStarHistograms ExecuteFitPhiWrongAsKStarShapes ExecuteFitKShortWrongAsKStarShapes

ExecuteMakeKStarSignalOnlyHistograms: MakeKStarSignalOnlyHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarSignalOnlyHistograms.cpp \
		-o ExecuteMakeKStarSignalOnlyHistograms \
		$(ROOTLIBS)

ExecuteFitKStarSignalOnlyShapes: FitKStarSignalOnlyShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSignalOnlyShapes.cpp \
		-o ExecuteFitKStarSignalOnlyShapes \
		$(ROOTLIBS)

ExecuteMakeKStarWrongAssignmentHistograms: MakeKStarWrongAssignmentHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarWrongAssignmentHistograms.cpp \
		-o ExecuteMakeKStarWrongAssignmentHistograms \
		$(ROOTLIBS)

ExecuteMakePhiWrongAsKStarHistograms: MakePhiWrongAsKStarHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakePhiWrongAsKStarHistograms.cpp \
		-o ExecuteMakePhiWrongAsKStarHistograms \
		$(ROOTLIBS)

ExecuteMakeKShortWrongAsKStarHistograms: MakeKShortWrongAsKStarHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKShortWrongAsKStarHistograms.cpp \
		-o ExecuteMakeKShortWrongAsKStarHistograms \
		$(ROOTLIBS)

ExecuteCombineKStarAssignments: CombineKStarAssignments.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		CombineKStarAssignments.cpp \
		-o ExecuteCombineKStarAssignments \
		$(ROOTLIBS)

ExecuteMakeKStarSBHistograms: MakeKStarSBHistograms.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		MakeKStarSBHistograms.cpp \
		-o ExecuteMakeKStarSBHistograms \
		$(ROOTLIBS)

ExecuteFitKStarSB: FitKStarSB.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSB.cpp \
		-o ExecuteFitKStarSB \
		$(ROOTLIBS)

ExecuteFitKStarSBWithCrossFeed: FitKStarSBWithCrossFeed.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSBWithCrossFeed.cpp \
		-o ExecuteFitKStarSBWithCrossFeed \
		$(ROOTLIBS)

ExecuteFitKStarSBFloatWidth: FitKStarSBFloatWidth.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKStarSBFloatWidth.cpp \
		-o ExecuteFitKStarSBFloatWidth \
		$(ROOTLIBS)

ExecuteFitPhiWrongAsKStarShapes: FitPhiWrongAsKStarShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitPhiWrongAsKStarShapes.cpp \
		-o ExecuteFitPhiWrongAsKStarShapes \
		$(ROOTLIBS)

ExecuteFitKShortWrongAsKStarShapes: FitKShortWrongAsKStarShapes.cpp
	g++ -O3 \
		$(ROOTCFLAGS) \
		FitKShortWrongAsKStarShapes.cpp \
		-o ExecuteFitKShortWrongAsKStarShapes \
		$(ROOTLIBS)
