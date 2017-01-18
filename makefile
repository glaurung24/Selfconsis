CC = g++
CFLAGS = -Wall -std=c++11 -fopenmp
CFLAGS_DEBUG = -Wall -std=c++11 -fopenmp -g
TBTK_DIR = ${TBTK_dir}
#CC = nvcc
#CFLAGS = -std=c++11 --compiler-options "-fopenmp -O3 -Wall"

all:
	$(CC) $(CFLAGS) src/main.cpp src/Calculation.cpp src/ProcessArgs.cpp -Iinclude -I$(TBTK_DIR)/hdf5/hdf5-build/include -L$(TBTK_DIR)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTBTK -lblas -llapack -lhdf5 -lhdf5_cpp -lcusparse -I$(TBTK_DIR)/hdf5/hdf5-build/hdf5/include/
#-I$(TBTK_DIR)/hdf5/hdf5-1.8.16/include/

clean:
	rm -r build/*

debug:
	$(CC) $(CFLAGS_DEBUG) src/main.cpp src/Calculation.cpp src/ProcessArgs.cpp -Iinclude -I$(TBTK_DIR)/hdf5/hdf5-build/include -L$(TBTK_DIR)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTBTK -lblas -llapack -lhdf5 -lhdf5_cpp -lcusparse -I$(TBTK_DIR)/hdf5/hdf5-build/hdf5/include/ -I$(TBTK_DIR)/hdf5/hdf5-1.8.16/include/


