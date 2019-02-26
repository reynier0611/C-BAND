ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../Hipo
HIPOLIBS    := -L../lib -lhipo

LZ4LIBS     := -llz4

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)


all: example

example: example.o
	$(CXX) -o example $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o readParticles readBanks analysis clas12event_example *~

%.o: %.cc
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS)
