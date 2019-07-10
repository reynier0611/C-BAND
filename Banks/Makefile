ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

CXX       := g++
CXXFLAGS  += -std=c++11 -Wall -fPIC $(ROOTCFLAGS)
LD        := g++
LDFLAGS   := $(ROOTLDFLAGS)

SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

LIBDIR  = ../lib

HIPOINCLUDES := -I../Hipo

libhipo: checkdirs $(OBJECTS)
	ar cr ../lib/libclas12banks.a $(OBJECTS)
#       $(CXX) -shared -fPIC  -o ../lib/libclas12banks.so $(OBJECTS) -llz4

show:
	@echo $(OBJECTS)

checkdirs: $(LIBDIR)

$(LIBDIR):
	@mkdir -p $(LIBDIR)

clean:
	@echo 'Removing the build files....'
	@rm -rf *~ *.o

%.o: %.cpp
	$(CXX) -fPIC -c $< -O2 -D__LZ4__ $(HIPOINCLUDES) $(ROOTCFLAGS)

