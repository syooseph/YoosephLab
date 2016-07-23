CPP=g++
CURRENT=/home/cczhong/Codes/GRASPxp_current
SUFSORT=$(CURRENT)/Libraries/libdivsufsort-2.0.1
BOOST_INC=/usr/local/include
BOOST_LIB=/usr/local/lib
THREADPOOL=$(CURRENT)/Libraries/threadpool
CC = $(CPP)
CFLAGS = -Wall -Wno-deprecated -std=c++0x
OPT = -O3
INCPATH = -I$(SUFSORT)/include -I$(THREADPOOL) -I$(BOOST_INC)
LIBPATH = -L$(SUFSORT)/lib -L$(BOOST_LIB)
LIBS = -Wl,-Bstatic $(SUFSORT)/lib/libdivsufsort64.a -Wl,-Bdynamic -lboost_iostreams -lboost_filesystem -lboost_system -lboost_thread -lboost_program_options -lrt -lpthread -lm -lz
LIBS_STATIC = -Wl,-Bstatic $(SUFSORT)/lib/libdivsufsort64.a /usr/local/lib/boost_1_55_0/stage/lib/libboost_iostreams.a /usr/local/lib/boost_1_55_0/stage/lib/libboost_filesystem.a /usr/local/lib/boost_1_55_0/stage/lib/libboost_system.a /usr/local/lib/boost_1_55_0/stage/lib/libboost_program_options.a /usr/local/lib/boost_1_55_0/stage/lib/libboost_thread.a -Wl,-Bdynamic  -lrt -lpthread -lm -lz
OBJS_ALL = main_build.o main_map.o main_assemble_p.o remap.o seq_align_extend_p.o hmm_profile.o contig_refinement_p.o assemble_extend_p.o reachable_reads.o database_index.o reduced_alphabet.o sequence_build.o gsa.o sfa.o interval_array.o
OBJS_XP = main_assemble_p.o seq_align_extend_p.o hmm_profile.o contig_refinement_p.o assemble_extend_p.o reachable_reads.o database_index.o reduced_alphabet.o sequence_build.o gsa.o sfa.o interval_array.o
OBJS_BUILD = main_build.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o
OBJS_MAP = main_map.o remap.o assemble_extend.o reachable_reads.o database_index.o sequence_build.o reduced_alphabet.o gsa.o sfa.o interval_array.o


all:	$(OBJS_ALL)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o graspxp-build $(LIBS_STATIC)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_XP) -o graspxp-assemble $(LIBS_STATIC)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_MAP) -o graspxp-map $(LIBS_STATIC)
graspxp-assemble: $(OBJS_XP)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_XP) -o graspxp-assemble $(LIBS_STATIC) 
graspxp-map: $(OBJS_MAP)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_MAP) -o graspxp-map $(LIBS_STATIC)
graspxp-build: $(OBJS_BUILD)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o graspxp-build $(LIBS_STATIC)
main_assemble_p.o: main_assemble_p.cc assemble_functor_p.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_assemble_p.cc
main_map.o: main_map.cc
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_map.cc
main_build.o: main_build.cc
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_build.cc
assemble_extend_p.o: assemble_extend_p.cc assemble_extend_p.h seq_align_extend_p.h hmm_profile.h reachable_reads.h database_index.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c assemble_extend_p.cc
remap.o: remap.cc remap.h assemble_extend.h gsa.h sequence_build.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c remap.cc
seq_align_extend_p.o: seq_align_extend_p.cc seq_align_extend_p.h seq_align.h hmm_profile.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c seq_align_extend_p.cc
contig_refinement_p.o: contig_refinement_p.cc contig_refinement_p.h sequence_build.h assemble_extend_p.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c contig_refinement_p.cc
hmm_profile.o: hmm_profile.cc hmm_profile.h reachable_reads.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c hmm_profile.cc
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
