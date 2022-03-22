CXX := g++
CXXFLAGS := -Wall -Wextra -std=c++11 -g

MKLPATH=/opt/intel/oneapi/mkl/latest/lib/intel64

INCLUDE := -I/usr/include/hdf5 -I/usr/include/hdf5/serial -I/usr/local/include/skynet
#LDLIBS := -lpardiso600-GNU720-X86-64
LDLIBS := -L${MKLPATH} -lmkl_core -lmkl_sequential -lmkl_intel_lp64\
-lSkyNet -lhdf5_cpp -lhdf5_serial -lgsl \
-lboost_system -lboost_filesystem -lboost_serialization -lblas -lgslcblas \
-llapack -fopenmp -lcairomm-1.0

all: bsm simple

bsm: src/r-process.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDLIBS) $^ -o $@ 

simple: src/r-canonical.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDLIBS) $^ -o $@ 

clean:
	rm -f bsm simple history.* final_abundance.txt 
