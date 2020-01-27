ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../../Hipo4 -I../../Banks4 -I../../Physics
HIPOLIBS    := -L../../lib -lhipo4 -lclas12banks4 -lclasphysics

LZ4LIBS     := -L/group/clas12/packages/lz4/lib -llz4
LZ4INCLUDES := -I/group/clas12/packages/lz4/lib

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

all: reactions reactions_BAND

reactions: reactions.o
	$(CXX) -o reactions $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

reactions_BAND: reactions_BAND.o
	$(CXX) -o reactions_BAND $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o reactions reactions_BAND *~

%.o: %.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
