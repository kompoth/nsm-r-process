CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++11

INCLUDE := -I/usr/include/hdf5 -I/usr/include/hdf5/serial -I/usr/local/include/skynet -Iext
LDLIBS := -lSkyNet -lhdf5_cpp -lpardiso600-GNU720-X86-64 -lhdf5_serial -lgsl \
-lboost_system -lboost_filesystem -lboost_serialization -lblas -lgslcblas \
-llapack -fopenmp -lcairomm-1.0

EXT_FILES := ext/*.cpp 

all: run test

run: r-process.cpp ${EXT_FILES} 
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDLIBS) $^ -o $@ 

test: test.cpp ${EXT_FILES}
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDLIBS) $^ -o $@ 

clean:
	rm -f run test
