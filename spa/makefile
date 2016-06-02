#C++      = $(CPP)
C++      = /usr/local/packages/gcc-4.7.1/bin/g++
CFLAGS   = -g -Wall -Wno-deprecated
OPT      = -O2
#INCPATH  = -I$(SEQAN)
INCPATH  = -I/usr/local/projdata/0599/projects/SPA/yyang/packages/Seqan_v1.3
LIBS     =  -lrt -lm -lz /usr/local/packages/boost-1.51.0/lib/libboost_iostreams.a  /usr/local/packages/boost-1.51.0/lib/libboost_program_options.a
OBJS0    = postfgs.o
OBJS1    = preprocess.o
OBJS2    = graph.o semiglobal.o galignment.o vprofile.o path.o msa.o assembly.o main.o
OBJS3    = rtran.o
EXE0	 = postfgs
EXE1     = prespa
EXE2     = spa
EXE3     = rtran
LONG     = 0
DEV      = 0

.cpp.o:
	$(C++) $(CFLAGS) -c $< $(INCPATH) -DLONGKMER=$(LONG) -DDEVELOPER=$(DEV)

install: postfgs prespa spa rtran

postfgs: $(OBJS0)
	$(C++)  $(OBJS0) $(OPT) -o $(EXE0) $(LIBS)
prespa: $(OBJS1)
	$(C++)  $(OBJS1) $(OPT) -o $(EXE1) $(LIBS)  
spa: $(OBJS2)
	$(C++) $(OBJS2) $(OPT) -o $(EXE2) $(LIBS)  
rtran: $(OBJS3)
	$(C++)  $(OBJS3) $(OPT) -o $(EXE3) $(LIBS)  

all: clean install

clean:
	rm -f *.o *core
