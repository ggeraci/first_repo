OBJS = simulation.o proc.o params.o gridsize.o particle.o tensor0.o tensor1.o  poisson.o poisson_IO.o communicator.o grid.o
CPP = mpicxx
CFLAGS = -O3
FFT_DIR=/home/ggeraci/TOOLS/fft
FFTW_DIR=/home/ggeraci/TOOLS/FFTW_2.1.5
HYPRE_DIR = /home/ggeraci/TOOLS/HYPRE
CINCLUDES = -I$(HYPRE_DIR)/include -I$(FFTW_DIR)/include
CDEFS     = -DHAVE_CONFIG_H -DHYPRE_TIMING
CFLAGS2 = -O3 -DFFT_FFTW  $(CINCLUDES) $(CDEFS)
LFLAGS = -L$(FFT_DIR)/Obj_certainty/ -lfft -L$(FFTW_DIR)/lib -lfftw -lfftw_mpi -lfftw -lm -L$(HYPRE_DIR)/lib -lHYPRE -lstdc++


default: box

box: $(OBJS)
	$(CPP) $(OBJS) -o box $(LFLAGS)

simulation.o: simulation.cpp proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp communicator.h communicator.cpp grid.h grid.cpp particle.h particle.cpp
	$(CPP) $(CFLAGS2) -c simulation.cpp

grid.o: proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp poisson.h poisson.cpp poisson_IO.h poisson_IO.cpp communicator.h communicator.cpp grid.h grid.cpp particle.h particle.cpp
	$(CPP) $(CFLAGS2) -c grid.cpp

communicator.o: proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp communicator.h communicator.cpp
	$(CPP) $(CFLAGS) -c communicator.cpp

poisson.o: $(FFT_DIR)/fft_3d.h $(FFT_DIR)/fft_3d.c proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp poisson.h poisson.cpp
	$(CPP) $(CFLAGS2) -c poisson.cpp

poisson_IO.o: $(FFT_DIR)/fft_2d.h $(FFT_DIR)/fft_2d.c proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp poisson_IO.h poisson_IO.cpp
	$(CPP) $(CFLAGS2) -c poisson_IO.cpp

particle.o: communicator.h communicator.cpp gridsize.h gridsize.cpp proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp particle.h particle.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp
	$(CPP) $(CFLAGS) -c particle.cpp

tensor1.o: communicator.h communicator.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp
	$(CPP) $(CFLAGS) -c tensor1.cpp

tensor0.o: communicator.h communicator.cpp gridsize.h gridsize.cpp tensor0.h tensor0.cpp tensor1.h tensor1.cpp
	$(CPP) $(CFLAGS) -c tensor0.cpp

gridsize.o: proc.h proc.cpp params.h params.cpp gridsize.h gridsize.cpp
	$(CPP) $(CFLAGS) -c gridsize.cpp

params.o: params.h params.cpp
	$(CPP) $(CFLAGS) -c params.cpp

proc.o: proc.h proc.cpp
	$(CPP) $(CFLAGS) -c proc.cpp

clean:
	rm -rf *.o box