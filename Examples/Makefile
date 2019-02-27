ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../Hipo -I../Banks
HIPOLIBS    := -L../lib -lhipo -lclas12banks

LZ4LIBS     := -llz4

CXX       := g++
CXXFLAGS  := $(HIPOCFLAGS) $(ROOTLDFLAGS) $(ROOTCFLAGS)
#COMPILER_FLAGS = $(HIPOCFLAGS) $(CXXFLAGS) $(ROOTCFLAGS) $(ROOTLDFLAGS)
LDFLAGS   := $(ROOTLDFLAGS)
LD	  := g++

SRC_FILES   := $(wildcard *.cxx)
EXE_FILES   := $(SRC_FILES:%.cxx=%)

.PHONY: all clean

all: $(EXE_FILES)

%: %.o
	$(LD) $(LDFLAGS) $(HIPOLIBS) $(LZ4LIBS) $(ROOTLIBS) $^ -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS)  -c $< 

clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)
