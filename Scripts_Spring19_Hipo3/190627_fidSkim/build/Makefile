ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../../../Hipo -I../../../Banks
HIPOLIBS    := -L../../../lib -lhipo -lclas12banks

LZ4LIBS     := -L/group/clas12/packages/lz4/lib -llz4
LZ4INCLUDES := -I/group/clas12/packages/lz4/lib

INCLUDES    := -I../include

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

all: pcal_fiducial

pcal_fiducial: pcal_fiducial.o
	$(CXX) -o pcal_fiducial $< $(ROOTCFLAGS) $(ROOTLDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS)

clean:
	@echo 'Removing all build files'
	@rm -rf *.o pcal_fiducial *~

%.o: ../%.cpp
	$(CXX) -c $< -O2 $(ROOTCFLAGS) $(HIPOCFLAGS) $(INCLUDES) $(LZ4INCLUDES)
