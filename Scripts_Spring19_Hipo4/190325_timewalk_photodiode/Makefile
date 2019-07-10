ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../../Hipo4 -I../../Banks4
HIPOLIBS    := -L../../lib -lhipo4 -lclas12banks4

LZ4LIBS     := -L/group/clas12/packages/lz4/lib -llz4
LZ4INCLUDES := -I/group/clas12/packages/lz4/lib

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

all: test_parameters timeWalk_photodiode_perPMT

test_parameters: test_parameters.o
	$(CXX) -o test_parameters $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

timeWalk_photodiode_perPMT: timeWalk_photodiode_perPMT.o
	$(CXX) -o timeWalk_photodiode_perPMT $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o test_parameters timeWalk_photodiode_perPMT *~

%.o: %.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
