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

all: parameter_extraction parameter_test mean_time

parameter_extraction: parameter_extraction.o
	$(CXX) -o parameter_extraction $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

parameter_test: parameter_test.o
	$(CXX) -o parameter_test $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

mean_time: mean_time.o
	$(CXX) -o mean_time $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o parameter_extraction parameter_test mean_time *~

%.o: %.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(LZ4INCLUDES)
