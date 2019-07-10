ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../../Hipo -I../../Banks
HIPOLIBS    := -L../../lib -lhipo -lclas12banks

LZ4LIBS     := -llz4

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

all: mean_time

mean_time: mean_time.o
	$(CXX) -o mean_time $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o mean_time *~

%.o: %.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS)
