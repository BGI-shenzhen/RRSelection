
CXX=g++
CXXFLAGS=  -g -O2
BIN := ../bin
LDFLAGS=-lz
INCLUDE=  -L/usr/lib64/   -L/usr/lib/ -L./include/zlib/
all: $(BIN)/RRSelection

$(BIN)/RRSelection: $(BIN)/RRSelection.o 
	$(CXX)   $^ -o $@   $(LDFLAGS)  $(INCLUDE) 

$(BIN)/%.o: %.cpp
	$(CXX)  -c $(CXXFLAGS) $< -o $@    $(INCLUDE) 

clean:
	$(RM) $(BIN)/*.o
