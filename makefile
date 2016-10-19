CC = g++
CFLAGS = -Wall -std=c++11 -fopenmp
CFLAGS_DEBUG = -Wall -std=c++11 -fopenmp -g
TBTK_DIR = /home/andi/Documents/MasterThesis/TBTK
#CC = nvcc
#CFLAGS = -std=c++11 --compiler-options "-fopenmp"

all:
	$(CC) $(CFLAGS) src/main.cpp src/Calculation.cpp -Iinclude -I$(TBTK_DIR)/TBTK/calc/TightBindingLib/include -I$(TBTK_DIR)/hdf5/hdf5-build/include -L$(TBTK_DIR)/TBTK/calc/TightBindingLib/build -L$(TBTK_DIR)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTightBinding -lblas -llapack -lhdf5 -lhdf5_cpp -lcusparse

clean:
	rm -r build/*

debug:
	$(CC) $(CFLAGS_DEBUG) src/main.cpp src/Calculation.cpp -Iinclude -I$(TBTK_DIR)/TBTK/calc/TightBindingLib/include -I$(TBTK_DIR)/hdf5/hdf5-build/include -L$(TBTK_DIR)/TBTK/calc/TightBindingLib/build -L$(TBTK_DIR)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTightBinding -lblas -llapack -lhdf5 -lhdf5_cpp -lcusparse


