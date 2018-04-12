CUDACC = nvcc
MPICXX=mpicxx
CXX=g++
CUDAFLAGS=
#CUDAFLAGS=--maxrregcount=128
DEBUGFLAGS=-g
CUDALDFLAGS=-L/usr/local/cuda/lib64/ -lcuda -lcudart -lm -lnvToolsExt
LDFLAGS=-fopenmp
FILES = params.o main.o read_particles.o init.o mpi_shortcut.o archAPI.o memory_control.o
OUT_EXE = all

all: $(FILES)
	$(MPICXX) -o $(OUT_EXE) $(FILES) $(LDFLAGS) 

main.o:
	$(MPICXX) -c $(CUDAFLAGS) $(DEBUGFLAGS) main.cpp  
	
read_particles.o:
	$(MPICXX) -c $(DEBUGFLAGS) read_particles.cxx


mpi_shortcut.o: 
	$(MPICXX) -c $(DEBUGFLAGS) mpi_shortcut.cxx
	
params.o:
	$(MPICXX) -c $(DEBUGFLAGS) params.cxx
	
archAPI.o:
	$(MPICXX) -c $(DEBUGFLAGS) archAPI.cxx	
	
memory_control.o:
	$(MPICXX) -c $(DEBUGFLAGS) memory_control.cxx	

init.o:
	$(CXX) -c $(DEBUGFLAGS) init.cpp        

clean:
	rm -f *.o core $(OUT_EXE)

rebuild: clean build
