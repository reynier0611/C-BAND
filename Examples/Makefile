ROOTCFLAGS  := $(shell root-config --cflags)
ROOTLDFLAGS := $(shell root-config --ldflags)
ROOTLIBS    := $(shell root-config --libs) -lEG
ROOTGLIBS   := $(shell root-config --glibs)

HIPOCFLAGS  := -I../Hipo -I../Banks
HIPOLIBS    := -L../lib -lhipo -lclas12banks

LZ4FLAGS    := -llz4

CXX         := g++
CXXFLAGS    := $(HIPOCFLAGS) $(ROOTLDFLAGS) $(ROOTCFLAGS)
LDFLAGS     := $(ROOTLDFLAGS)
LD	    := g++

SRC_FILES   := $(wildcard *.cxx)
EXE_FILES   := $(SRC_FILES:%.cxx=%)

.PHONY: all clean

all: $(EXE_FILES)

%: %.o
	$(LD) $(LDFLAGS) $(HIPOLIBS) $(LZ4FLAGS) $(ROOTLIBS) $^ -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS)  -c $< 

clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)
