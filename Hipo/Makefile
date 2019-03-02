
SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)

LIBDIR  = ../lib

LZ4INCLUDES := -I/group/clas12/packages/lz4/lib/

CXXFLAGS = -std=c++11

libhipo: checkdirs $(OBJECTS)
	ar cr ../lib/libhipo.a $(OBJECTS)
#	$(CXX) -shared -fPIC  -o ../lib/libhipo.so $(OBJECTS) -llz4


show:
	@echo $(OBJECTS)

checkdirs: $(LIBDIR)

$(LIBDIR):
	@mkdir -p $(LIBDIR)

clean:
	@echo 'Removing the build files....'
	@rm -rf *~ *.o

%.o: %.cpp
	$(CXX) -fPIC -c $< -O2 -D__LZ4__ $(CXXFLAGS) $(LZ4INCLUDES)
