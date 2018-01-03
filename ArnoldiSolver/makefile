CC = g++
CFLAGS = -Wall -std=c++11 -fopenmp
#CC = nvcc
#CFLAGS = -std=c++11 --compiler-options "-fopenmp"

all:
	@echo "Building: BasicArnoldi"
	@$(CC) $(CFLAGS) src/main.cpp src/ProcessArgs.cpp src/Calculation.cpp -Iinclude -I$(TBTK_dir)/hdf5/hdf5-build/include -L$(TBTK_dir)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTBTK -lblas -llapack -lhdf5 -lhdf5_cpp -lsuperlu -larpack -lblas -I/usr/include/superlu -I$(TBTK_dir)/hdf5/hdf5-build/hdf5/include

clean:
	rm -r build/*

