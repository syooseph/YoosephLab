C++      = g++
CFLAGS   = -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT      = -O2

SPAHOME  = $(SPA_HOME)
SEQAN    = $(SPAHOME)/3rd/seqan
SUFSORT  = $(SPAHOME)/3rd/libdivsufsort-2.0.1/build
INCPATH  = -I$(SEQAN) -I$(SUFSORT)/include -I$(BOOST)/include
LIBPATH  = -L$(BOOST)/lib -L$(SUFSORT)/lib -L/usr/local/lib -L/usr/local/lib64
LIBS     = -lrt -lm -lz -lgomp -lpthread -lboost_iostreams -lboost_program_options -ldivsufsort

OBJS0    = param.o cmdargs_simple.o codon.o postfgs.o
OBJS1    = param.o cmdargs_simple.o part.o
OBJS2    = sfa.o gsa.o IntervalArray.o param.o cmdargs_simple.o preprocess.o
OBJS3    = meminfo.o smath.o graph.o IntervalArray.o sfa.o gsa.o ReadStartCount.o path.o semiglobal.o alignsummary.o PathKmerIndex.o PathKmerPosIndex.o KmerPathIndex.o PathEntry.o param.o loader.o extracter.o sectioner.o msectioner.o esectioner.o cluster.o merger.o  extender.o ConnectEntries.o connecter.o ljoiner.o sjoiner.o rjoiner.o rjoinerSA.o seed.o rcmat.o rplace.o ralign.o palign.o msa.o profile.o anchor.o metapath.o placer.o recruiter.o reporter.o assembler.o cmdargs_simple.o log.o spamain.o
OBJS4    = param.o cmdargs_simple.o semiglobal.o alignsummary.o path.o rcmat.o ReadStartCount.o ralign.o palign.o msa.o profile.o KmerPathIndex.o PathEntry.o sectioner.o anchor.o metapath.o placer.o reverser.o extender.o IntervalArray.o sfa.o gsa.o codon.o rtran.o log.o rplace.o

EXE0	 = $(SPAHOME)/bin/postfgs
EXE1     = $(SPAHOME)/bin/part
EXE2     = $(SPAHOME)/bin/prespa
EXE3     = $(SPAHOME)/bin/spa
EXE4     = $(SPAHOME)/bin/rtran

.cpp.o:
	$(C++) $(CFLAGS) $(INCPATH) -c $<

install:  postfgs part prespa spa rtran 

postfgs: $(OBJS0)
	$(C++) $(OPT) -o $(EXE0) $(OBJS0) $(LIBPATH) $(LIBS) 
	@echo
part: $(OBJS1)
	$(C++) $(OPT) -o $(EXE1) $(OBJS1) $(LIBPATH) $(LIBS) 
	@echo
prespa: $(OBJS2)
	$(C++) $(OPT) -o $(EXE2) $(OBJS2) $(LIBPATH) $(LIBS) 
	@echo
spa: $(OBJS3)
	$(C++) $(OPT) -o $(EXE3) $(OBJS3) $(LIBPATH) $(LIBS) 
	@echo
rtran: $(OBJS4)
	$(C++) $(OPT) -o $(EXE4) $(OBJS4) $(LIBPATH) $(LIBS) 
	@echo

all: clean install

clean:
	rm -f $(SPAHOME)/bin/* *.o *core

