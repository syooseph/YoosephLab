CPP=g++
CURRENT=/home/cczhong/Codes/Grasp
SUFSORT=$(CURRENT)/Libraries/libdivsufsort-2.0.1
BOOST_INC=/usr/local/include
BOOST_LIB=/usr/local/lib
THREADPOOL=$(CURRENT)/Libraries/threadpool
CC = $(CPP)
CFLAGS = -Wall -Wno-deprecated -std=c++0x
OPT = -O3
INCPATH = -I$(SUFSORT)/include -I$(THREADPOOL) -I$(BOOST_INC)
LIBPATH = -L$(SUFSORT)/lib -L$(BOOST_LIB)
LIBS = -lboost_iostreams -lboost_filesystem -lboost_system -lboost_thread -lboost_program_options -ldivsufsort_la-divsufsort -ldivsufsort_la-sssort -ldivsufsort_la-trsort -ldivsufsort_la-utils -lrt -lpthread -lm -lz
OBJS_ALL = main_build.o main_search.o guided_assemble.o extend_and_assemble.o assembly_graph.o index_sample.o sfa.o gsa.o interval_array.o
OBJS_BUILD = main_build.o guided_assemble.o extend_and_assemble.o assembly_graph.o index_sample.o sfa.o gsa.o interval_array.o
OBJS_SEARCH = main_search.o guided_assemble.o extend_and_assemble.o assembly_graph.o index_sample.o sfa.o gsa.o interval_array.o
All:	$(OBJS_ALL)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o Build $(LIBS)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_SEARCH) -o Search $(LIBS)
Build:	$(OBJS_BUILD)
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_BUILD) -o Build $(LIBS)
Search: $(OBJS_SEARCH) 
	$(CC) $(INCPATH) $(LIBPATH) $(CFLAGS) $(OPT) $(OBJS_SEARCH) -o Search $(LIBS)
main_build.o: main_build.cc sequence.h extend_and_assemble.h index_sample.h gsa.h guided_assemble.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_build.cc
main_search.o: main_search.cc sequence.h extend_and_assemble.h index_sample.h gsa.h guided_assemble.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c main_search.cc
guided_assemble.o: guided_assemble.cc extend_and_assemble.h extend_functor.h certificate_functor.h realignment_functor.h guided_assemble.h index_sample.h tree_basic.h gsa.h sequence.h timer.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c guided_assemble.cc
extend_and_assemble.o: extend_and_assemble.cc extend_and_assemble.h tree_basic.h seq_align.h seq_align_extend.h scoring_function.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c extend_and_assemble.cc
assembly_graph.o: assembly_graph.cc assembly_graph.h index_sample.h sfa.h gsa.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c assembly_graph.cc
index_sample.o: index_sample.cc index_sample.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c index_sample.cc
gsa.o: gsa.cc gsa.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c gsa.cc
sfa.o: sfa.cc sfa.h timer.h file.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c sfa.cc
interval_array.o: interval_array.cc interval_array.h
	$(CC) $(INCPATH) $(CFLAGS) $(OPT) -c interval_array.cc

new:
	-touch *.cc
	-make Build
	-make Search
clean:
	-rm Build Search $(OBJS_ALL)
