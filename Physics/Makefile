ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

#PParticle.o: PParticle.cpp PParticle.h
#	$(CXX) PParticle.cpp -c -o PParticle.o $(ROOTCFLAGS)

#OBJECTS := PParticle.o

LIBDIR  = ../lib

libPhysics: checkdirs $(OBJECTS)
	ar cr ../lib/libclasphysics.a $(OBJECTS)

show:
	@echo $(OBJECTS)

checkdirs: $(LIBDIR)

$(LIBDIR):
	@mkdir -p $(LIBDIR)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o

%.o: %.cpp
	$(CXX) -fPIC -c $< -O2 $(ROOTCFLAGS) 
