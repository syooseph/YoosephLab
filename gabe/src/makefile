## To compile,
## $ make
## By default, it uses a sparse matrix. To use full matrix,
## $ make SPARSE=0
## 

C++      = g++
CFLAGS   = -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT      = -O2

LIBS     = -lrt -lz -lgomp 
#LIBS     = -lrt -lz -static-libstdc++ -static-libgcc /usr/local/packages/gcc-4.9.3/lib64/libgomp.a

OBJS0    = bwa_reader.o genome_matrix.o preprocessor.o
OBJS1    = matrix.o mg.o gem.o gem_handler.o

EXE0     = ../bin/gabe_pre
EXE1     = ../bin/gabe_run

SPARSE   = 1

.PHONY: all clean install

.cpp.o:
	$(C++) $(CFLAGS) -c $< -DSPARSE=$(SPARSE) 

install:  pre gem

pre: $(OBJS0)
	$(C++) $(OPT) -o $(EXE0) $(OBJS0) $(LIBPATH) $(LIBS)	
gem: $(OBJS1)
	$(C++) $(OPT) -o $(EXE1) $(OBJS1) $(LIBPATH) $(LIBS) 

all: clean install

clean:
	rm -f $(EXE0) $(EXE1) *.o *core
