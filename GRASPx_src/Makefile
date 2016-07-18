# Please set the following paths before compilation

# This should be the path to the home directory of GRASPx
CURRENT=/home/cczhong/Codes/GRASPx_current
# This should be the home directory of BOOST library; remember also adding lib directory into LD_LIBDARY_PATH
BOOST=/usr/local

# Inferring the paths
CPP=g++
SUFSORT=$(CURRENT)/Libraries/libdivsufsort-2.0.1
BOOST_INC=$(BOOST)/include
BOOST_LIB=$(BOOST)/lib
THREADPOOL=$(CURRENT)/Libraries/threadpool
CC = $(CPP)
CFLAGS = -Wall -Wno-deprecated -std=c++0x
OPT = -O3
INCPATH = -I$(SUFSORT)/include -I$(THREADPOOL) -I$(BOOST_INC)
LIBPATH = -L$(SUFSORT)/lib -L$(BOOST_LIB)
LIBS = -Wl,-Bstatic $(SUFSORT)/lib/libdivsufsort64.a -Wl,-Bdynamic -lboost_iostreams -lboost_filesystem -lboost_system -lboost_thread -lboost_program_options -lrt -lpthread -lm -lz
OBJS_ALL = main_build.o main_assemble.o main_map.o remap.o contig_refinement.o greedy_assembly.o assemble_extend.o read_alignment.o reachable_reads.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_BUILD = main_build.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_ASSEMBLE = main_assemble.o contig_refinement.o assemble_extend.o read_alignment.o reachable_reads.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_MAP = main_map.o remap.o assemble_extend.o reachable_reads.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_TEMP = main_all.o remap.o contig_refinement.o greedy_assembly.o assemble_extend.o read_alignment.o reachable_reads.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_UNITIGER = unitiger.o database_index.o sequence_build.o main_unitiger.o gsa.o sfa.o interval_array.o reduced_alphabet.o
all:	$(OBJS_ALL)
#	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_TEMP) -o Build $(LIBS)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o graspx-build $(LIBS)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_ASSEMBLE) -o graspx-assemble $(LIBS)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_MAP) -o graspx-map $(LIBS)
unitiger: $(OBJS_UNITIGER)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_UNITIGER) -o unitiger $(LIBS) 
graspx-map: $(OBJS_MAP)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_MAP) -o graspx-map $(LIBS)
graspx-assemble: $(OBJS_ASSEMBLE)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_ASSEMBLE) -o graspx-assemble $(LIBS)
graspx-build: $(OBJS_BUILD)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o graspx-build $(LIBS)
#main_all.o: main_all.cc
#	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_all.cc
main_map.o: main_map.cc
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_map.cc
main_assemble.o: main_assemble.cc assemble_functor.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_assemble.cc
main_build.o: main_build.cc
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_build.cc
main_unitiger.o: main_unitiger.cc
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_unitiger.cc
unitiger.o: unitiger.cc database_index.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c unitiger.cc
remap.o: remap.cc remap.h assemble_extend.h gsa.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c remap.cc
contig_refinement.o: contig_refinement.cc contig_refinement.h assemble_extend.h sequence_build.h seq_align.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c contig_refinement.cc
assemble_extend.o: assemble_extend.cc assemble_extend.h reachable_reads.h sequence_build.h scoring_function.h seq_align.h seq_align_extend.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c assemble_extend.cc
greedy_assembly.o:  greedy_assembly.cc greedy_assembly.h reachable_reads.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c greedy_assembly.cc
read_alignment.o: read_alignment.cc read_alignment.h reachable_reads.h scoring_function.h seq_align.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c read_alignment.cc
reachable_reads.o: reachable_reads.cc reachable_reads.h scoring_function.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c reachable_reads.cc
database_index.o: database_index.cc database_index.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c database_index.cc
reduced_alphabet.o: reduced_alphabet.cc reduced_alphabet.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c reduced_alphabet.cc
sequence_build.o: sequence_build.cc sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c sequence_build.cc 
gsa.o: gsa.cc gsa.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c gsa.cc
sfa.o: sfa.cc sfa.h timer.h file.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c sfa.cc
interval_array.o: interval_array.cc interval_array.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c interval_array.cc

new:
	-touch *.cc
	-make
clean:
	-rm $(OBJS_ALL)
