CC = g++
CFLAGS = -Wall -std=c++11 -fopenmp
INIT = /home/andi/Documents/MasterThesis/Lunarc/TBTK/init_session.sh
#CC = nvcc
#CFLAGS = -std=c++11 --compiler-options "-fopenmp"

all:
	$(INIT)
	$(CC) $(CFLAGS) src/main.cpp -I$(TBTK_dir)/TBTK/calc/TightBindingLib/include -I$(TBTK_dir)/hdf5/hdf5-build/include -L$(TBTK_dir)/TBTK/calc/TightBindingLib/build -L$(TBTK_dir)/hdf5/hdf5-build/hdf5/lib -o build/a.out -lTightBinding -lblas -llapack -lhdf5 -lhdf5_cpp -lcusparse

clean:
	rm -r build/*

