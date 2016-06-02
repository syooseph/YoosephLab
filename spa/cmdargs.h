/** 
 * \file cmdargs.h
 * Command line argument handler.
 * This module parses command line arguments.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \pre		  boost library must be installed. Also, compile with -lboost_program_options.
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __CMDARGS_H__
#define __CMDARGS_H__

#include <iostream>
#include <string>
#include <fstream>
#include <iterator>
#include <boost/program_options.hpp>

#include "default.h"
#include "param.h"

namespace po = boost::program_options;


/** 
 * Command line argument parser.
 */
namespace args 
{
	/**
	 * Usage message of post processing FragGeneScan program.
	 */
	const std::string USAGE_POSTFGS = 
		"Usage: postfgs [options]\n"
		" e.g.: postfgs -a fgs.faa -n fgs.ffn -o fgs.out reads1.fasta [reads2.fasta ...]\n";
	
	/**
	 * Usage message of preprocessing of sequences.
	 */
	const std::string USAGE_PREPARE = 
		"Usage: prerun [options]\n"
		" e.g.: prerun -k 6 -G graph -I index -i input.faa [input2.faa ...]\n";
	/**
	 * Usage message of main assembler
	 */
	const std::string USAGE_ASSEM = 
		"Usage: spa [options] inputs\n"
		" e.g.: spa -k 6 -G graph -I index --profile --alignment -M -i input.faa [input2.faa ...]\n";
	
	/**
	 * Usage message of post processing
	 */
	const std::string USAGE_POST_SPA = 
		"Usage: rtran [options] inputs\n"
		" e.g.: rtran -a fgs.faa -d fgs.ffn -p spa.place\n";

	const std::string USAGE_METISIN =
		"Usage: metisin options";

	const std::string USAGE_DEVIDE =
		"Usage: devide options";

	const std::string USAGE_MERGE =
		"Usage: merge options";

	/**
	 * Check whether valid k-mer size is given.
	 */
	void checkKmerSize(int k)
	{
#if LONGKMER == 1
		if ( k > 12 ) {
			std::cerr << "[Error] Maximum k-mer size is 12\n";
			exit(1);
		} 
#else
		if ( k > 6 ) {
			std::cerr << "[Error] Maximum k-mer size of the current binary is 6\n";
			std::cerr << "For the bigger kmer size, recompile with longer kmer support\n";
			exit(1);
		}
#endif
	}

	/**
	 * \param	filename	a file with a list of sequence files.
	 * \param   vector		empty vector.
	 * \post	Vector of sequence file names.
	 */
	void loadFileList(const char *listfile, std::vector<std::string> &input_files) 
	{
		std::ifstream inf(listfile, std::ifstream::in );

		if ( !inf ) {
			std::cerr << "[Error] Can't open " << listfile << "\n";
			exit(1);
		}

		std::string file;

		inf >> file;	
		while ( !inf.eof() ) {
			input_files.push_back( file );
			inf >> file;
		}
		inf.close();
	}

	/**
	 * Echo the command that a user typed.
	 */
	void printCommand(int ac, char **av) 
	{
		for ( int i = 0; i < ac; i++ )
			std::cout << av[i] << " ";
		std::cout << "\n";
	}

	/**
	 * Parse command line arguments for FragGeneScan post processing
	 */ 
    void parse_cmd_args_postfgs(int ac, char **av, 
								std::vector<std::string> &input_files, 
								std::string &fgs_out,
								std::string &fgs_ffn,
								std::string &fgs_faa,
								std::string &out_dir,
								bool &verbose)
	{
		try {
			std::string file_list;

			po::options_description generic("Generic options");
			generic.add_options()
				("help,h", "produce help message")
				;

			po::options_description list("File options");
			list.add_options()
				("out-dir,t", po::value<std::string>(&out_dir)->default_value(OUTDIR), "output directory")
				("input-list,l", po::value<std::string>(&file_list), "input file list")
				("input-file,i", po::value< std::vector<std::string> >(), "input file")
				("out-file,o", po::value<std::string>(&fgs_out), "FragGeneScan out file")
				("ffn-file,n", po::value<std::string>(&fgs_ffn), "FragGeneScan ffn file")
				("faa-file,a", po::value<std::string>(&fgs_faa), "FragGeneScan faa file")
				//("compressed,Z", "compressed")
				("verbose,V", "verbose switch on")
				;
		
			// Hidden options, will be allowed both on command line and
			// in config file, but will not be shown to the user.
			po::options_description hidden("Hidden options");
			hidden.add_options()
				("input-file(s)", po::value< std::vector<std::string> >(), "input file")
				;
		
			po::options_description cmdline_options;
			cmdline_options.add(generic).add(list).add(hidden);

			po::options_description visible("Options");
			visible.add(generic).add(list);
        
			po::positional_options_description p;
			p.add("input-file", -1);
        
			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).positional(p).run(), vm);
			notify(vm);


			if (vm.count("verbose")) 
				verbose = true;

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_POSTFGS << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if ( ! vm.count("out-file") ||
				 ! vm.count("ffn-file") ||
				 ! vm.count("faa-file") ) {
				std::cerr << "[Error] FGS files must be provided\n";
				std::cout << "\n" << USAGE_POSTFGS << "\n\n";
				std::cout << visible;
				exit(1);
			}

			if ( ! vm.count("input-list") && 
				 ! vm.count("input-file") ) {
				std::cerr << "[Error] Input reads must be given\n";
				std::cout << "\n" << USAGE_POSTFGS << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
			if ( vm.count("input-list")) {
				file_list = vm["input-list"].as<std::string>();
				loadFileList(file_list.c_str(), input_files);
			}
			if (vm.count("input-file")) {
				input_files = vm["input-file"].as< std::vector<std::string> >();
			}

		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    

	}

	/**
	 * Parse command line arguments for sequence pre-processing.
	 */ 
	void parse_cmd_args_prepare(int ac, 
								char **av, 
								int &k, 
								std::vector<std::string> &input_files, 
								std::string &graph_input,
								std::string &index_file,
								std::string &outdir )
	{
		try {
			std::string file_list;

			/* Declare a group of options that will be allowed only on command line */
			//po::options_description generic("Generic options");
			po::options_description generic("");
			generic.add_options()
				//("version,v", "print version string")
				("help,h", "This help message")
				//;

			//po::options_description parameter("Program parameters");
			//po::options_description parameter("");
			//parameter.add_options()
				("kmer,k", po::value<int>(&k)->default_value(KMER_SIZE), "Size of k-mer")
				//;

			//po::options_description list("File options");
			//po::options_description list("");
			//list.add_options()
				("input-file,i", po::value< std::vector<std::string> >(), "Input file(s)")
				("input-list,l", po::value<std::string>(&file_list), "File of multiple sequence locations")
				("graph-file,G", po::value<std::string>(&graph_input), "Graph file")
				("index-file,I", po::value<std::string>(&index_file), "Inverted index file")
				("out-dir,o", po::value<std::string>(&outdir)->default_value(OUTDIR), "Output directory")
				//;
				
				;
			
			/* Hidden options, will be allowed both on command line and
			 * in config file, but will not be shown to the user. */
			po::options_description hidden("Hidden options");
			hidden.add_options()
				("input-file(s)", po::value< std::vector<std::string> >(), "input file")
				;
		
			po::options_description cmdline_options;
			//cmdline_options.add(generic).add(parameter).add(list).add(hidden);
			cmdline_options.add(generic).add(hidden);

			po::options_description visible("Options");
			//visible.add(generic).add(parameter).add(list);
			visible.add(generic);
        
			po::positional_options_description p;
			p.add("input-file", -1);
        
			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).positional(p).run(), vm);
			notify(vm);

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_PREPARE << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("input-list") && 
				 ! vm.count("input-file") ) {
				std::cerr << "[Error] Input file(s) must be given\n";
				std::cout << "\n" << USAGE_PREPARE << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
			if ( vm.count("input-list")) {
				file_list = vm["input-list"].as<std::string>();
				loadFileList(file_list.c_str(), input_files);
			}
			if (vm.count("input-file")) {
				input_files = vm["input-file"].as< std::vector<std::string> >();
			}

			if ( ! vm.count("graph-file") ) {
				graph_input = outdir + "/gin.bin";
			}
			if ( ! vm.count("index-file") ) {
				index_file = outdir + "/ind.bin";
			}

			/* Error check */
			checkKmerSize(k);
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    
	}

/* 	/\** */
/* 	 * Parse command line arguements for assembler. */
/* 	 *\/  */
/* 	void parse_cmd_args_assem(int ac,  */
/* 							  char **av,  */
/* 							  int &k,  */
/* 							  int &min_depth, */
/* 							  int &min_share, */
/* 							  int &min_seed, */
/* 							  int &med_depth, */
/* 							  int &min_length, */
/* 							  int &back_trace, */
/* 							  int &base_depth, */
/* 							  int &gap_open, */
/* 							  int &gap_ext, */
/* 							  int &latch_length, */
/* 							  int &latch_support, */
/* 							  int &pairend_overlap, */
/* 							  int &overlap_support, */
/* 							  int &path_spur, */
/* 							  int &read_spur, */
/* 							  int &platform, */
/* 							  int &insert_size, */
/* 							  double &insert_sd, */
/* 							  double &read_align_score, */
/* 							  double &read_align_ratio, */
/* 							  double &filter_score, */
/* 							  double &merge_score, */
/* 							  double &latch_score, */
/* 							  double &latch_ratio, */
/* 							  double &kmer_percentile, */
/* 							  unsigned &seed, */
/* 							  std::vector<std::string> &input_files,  */
/* 							  std::string &graph_file,  */
/* 							  std::string &index_file, */
/* 							  std::string &out_dir, */
/* 							  std::string &file_suffix, */
/* 							  std::string &dump_suffix, */
/* 							  bool &pair_flag, */
/* 							  bool &trim_flag, */
/* 							  bool &path_flag, */
/* 							  bool &dump_flag, */
/* 							  bool &recruit_flag, */
/* 							  bool &finish_flag, */
/* 							  bool &latch_read_flag, */
/* 							  std::string &debug_file, */
/* 							  bool &latch_flag, */
/* 							  bool &profile_flag, */
/* 							  bool &align_flag, */
/* 							  bool &verbose) */
/* 	{ */
/* 		try { */
/* 			std::string file_list; */

/* 			po::options_description generic("Generic options"); */
/* 			generic.add_options() */
/* 				("version,v", "Print current version.") */
/* 				("help,h", "Produce help message.") */
/* 				; */

/* 			po::options_description graph("Graph options"); */
/* 			graph.add_options() */
/* 				("kmer-size,k", po::value<int>(&k)->default_value(KMER_SIZE), "Size of k-mer (node)") */
/* 				("min-coverage,c", po::value<int>(&min_depth)->default_value(MIN_DEPTH), "Minimum kmer coverage for graph trimming") */
/* 				("no-trim", "Skip trimming graph") */
/* 				; */
			
/* 			po::options_description path("Path search options"); */
/* 			path.add_options() */
/* 				("percentile,p", po::value<double>(&kmer_percentile)->default_value(KMER_PERCENTILE), "Seed k-mer: top x% percentile of kmers") */
/* 				("min-seed,s", po::value<int>(&min_seed)->default_value(SEED_START), "Minimum coverage of seed k-mer") */
/* 				("distance,d", po::value<int>(&back_trace)->default_value(BACK_TRACE), "Distance in path to find maximum read sharing node") */
/* 				("min-share,r", po::value<int>(&min_share)->default_value(SHARED_READS), "Minimum shared reads with neighboring node in path search") */
/* 				("min-length,n", po::value<int>(&min_length)->default_value(MIN_LENGTH), "Minimum path length") */
/* 				("no-path", "Skip path search") */
/* 				; */
			
/* 			po::options_description merge("Path merging/extension options"); */
/* 			merge.add_options() */
/* 				("long-overlap,L", po::value<int>(&latch_length)->default_value(LATCH_LENGTH), "Minimum overlapping length to extend a sequence pair") */
/* 				("long-support,S", po::value<int>(&latch_support)->default_value(LATCH_SUPPORT), "MInimum number of supporting reads in ovelapping paths") */
/* 				("merge-score,m", po::value<double>(&merge_score)->default_value(MERGE_SCORE), "Score for merging sequences") */
/* 				("extend-score,e", po::value<double>(&latch_score)->default_value(LATCH_SCORE), "Score for extending sequences") */
/* 				("path-spur,u", po::value<int>(&path_spur)->default_value(PATH_MISMATCH_ALLOW), "Maximum length of spur (path to path)") */
/* 				("no-merge", "Skip merging/extending paths") */
/* 				; */
			
/* 			po::options_description read("Read recruitment options"); */
/* 			read.add_options() */
/* 				("read-align,a", po::value<double>(&read_align_score)->default_value(READ_ALIGN_SCORE), "Minimum alignment score to be a part of a path") */
/* 				("read-ratio,y", po::value<double>(&read_align_ratio)->default_value(READ_ALIGN_RATIO), "Minimum length ratio of a read to be a part of a path") */
/* 				("read-spur,U", po::value<int>(&read_spur)->default_value(READ_MISMATCH_ALLOW), "Maximum length of spur (path to read)") */
/* 				("extend-reads", "Extend reads to a path") */
/* 				("no-recruit", "Skip read recruitment") */
/* 				; */

/* 			po::options_description pairend("Short path extension options"); */
/* 			pairend.add_options() */
/* 				("short-overlap,t", po::value<int>(&pairend_overlap)->default_value(PAIREND_OVERLAP), "Minimum overlapping length of short path extension") */
/* 				("short-support,z", po::value<int>(&overlap_support)->default_value(OVERLAP_SUPPORT), "Minimum supporting paired end reads for short overlap extension") */
/* 				("extend-ratio,E", po::value<double>(&latch_ratio)->default_value(LATCH_RATIO), "Minimum percentage of read length to be extended to a path") */
/* 				("no-short", "Skip short path extension") */
/* 				; */

/* 			po::options_description other("Other program options"); */
/* 			other.add_options() */
/* 				("med-coverage,C", po::value<int>(&med_depth)->default_value(MED_DEPTH), "Minimum median read coverage of amino acid in consensus") */
/* 				("base-depth,b", po::value<int>(&base_depth)->default_value(BASE_DEPTH), "Minimum base coverage of amino acid in a path") */
/* 				("filter-score,f", po::value<double>(&filter_score)->default_value(FILTER_SCORE), "Score for kmer filtering") */
/* 				("gap-open,g", po::value<int>(&gap_open)->default_value(GAP_OPEN), "Gap open penalty") */
/* 				("gap-ext,x", po::value<int>(&gap_ext)->default_value(GAP_EXT), "Gap extension penalty") */
/* 				("seed,R", po::value<unsigned>(&seed)->default_value(RAND_SEED), "random seed") */
/* 				("platform,T", po::value<int>(&platform)->default_value(PLATFORM), "Sequencing technology") */
/* 				("insert-size,w", po::value<int>(&insert_size)->default_value(INSERT_SIZE), "Inner insert size") */
/* 				("insert-sd,W", po::value<double>(&insert_sd)->default_value(INSERT_SD), "Insert std. deviation") */
/* 				("verbose,V", "verbose switch on") */
/* 				; */

/* 			po::options_description list("File options"); */
/* 			list.add_options() */
/* 				("input-file,i", po::value< std::vector<std::string> >(), "A sequence file") */
/* 				("input-list,l", po::value<std::string>(&file_list), "A list of sequence files") */
/* 				("graph-file,G", po::value<std::string>(&graph_file), "Graph file") */
/* 				("index-file,I", po::value<std::string>(&index_file), "Inverted index file") */
/* 				("output-dir,o", po::value<std::string>(&out_dir)->default_value(OUTDIR), "Output directory") */
/* 				("file-suffix,F", po::value<std::string>(&file_suffix), "File suffix for loading pre-written dump") */
/* 				("dump-suffix,X", po::value<std::string>(&dump_suffix), "File suffix for writing dump") */
/* 				("debug-file,B", po::value<std::string>(&debug_file), "K-mer list for debugging purpose") */
/* 				("dump-flag,D", "Dump objects") */
/* 				("pair-end,M", "Paired end reads") */
/* 				("profile,P", "Generate a profile file") */
/* 				("alignment,A", "Generate an alignment file") */
/* 				; */
		
/* 			// Hidden options, will be allowed both on command line and */
/* 			// in config file, but will not be shown to the user. */
/* 			po::options_description hidden("Hidden options"); */
/* 			hidden.add_options() */
/* 				("input-file(s)", po::value< std::vector<std::string> >(), "input file") */
/* 				; */

		
/* 			po::options_description cmdline_options; */
/* 			//cmdline_options.add(generic).add(parameter).add(list).add(hidden); */
/* 			cmdline_options.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list).add(hidden); */

/* 			po::options_description visible("Options"); */
/* 			//visible.add(generic).add(parameter).add(list); */
/* 			visible.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list); */
        
/* 			po::positional_options_description p; */
/* 			p.add("input-file", -1); */
        
/* 			po::variables_map vm; */
/* 			store(po::command_line_parser(ac, av). */
/* 				  options(cmdline_options).positional(p).run(), vm); */
/* 			notify(vm); */

/* 			if (vm.count("version")) { */
/* 				std::cerr << VERSION << "\n"; */
/* 				exit (1); */
/* 			} */
/* 			if (vm.count("help")) { */
/* 				std::cerr << "\n" << USAGE_ASSEM << "\n\n"; */
/* 				std::cout << visible << "\n"; */
/* 				exit (1); */
/* 			} */

/* 			if (vm.count("version")) { */
/* 				std::cout << "Multiple sources example, version 1.0\n"; */
/* 				exit (1); */
/* 			} */

			
/* 			if ( !vm.count("file-suffix") && !vm.count("graph-file") ) { */
/* 				std::cerr << "\nk-mer to read mapping file must be provided by either directly or dump\n"; */
/* 				std::cerr << "\n" << USAGE_ASSEM << "\n\n"; */
/* 				std::cerr << visible; */
/* 				exit(1); */
/* 			} */
/* 			if ( !vm.count("file-suffix") && !vm.count("index-file") ) { */
/* 				std::cerr << "\nk-mer to read mapping file must be provided by either directly or dump\n"; */
/* 				std::cerr << "\n" << USAGE_ASSEM << "\n\n"; */
/* 				std::cerr << visible; */
/* 				exit(1); */
/* 			} */

/* 			if (vm.count("pair-end")) */
/* 				pair_flag = true; */
/* 			if (vm.count("no-trim"))  */
/* 				trim_flag = false; */
/* 			if (vm.count("no-path"))  */
/* 				path_flag = false; */
/* 			if (vm.count("no-merge"))  */
/* 				latch_flag = false; */
/* 			if (vm.count("no-recruit"))  */
/* 				recruit_flag = false; */
/* 			if (vm.count("no-short"))  */
/* 				finish_flag = false; */
/* 			if (vm.count("extend-reads"))  */
/* 				latch_read_flag = true; */
/* 			if (vm.count("dump-flag"))  */
/* 				dump_flag = true; */
/* 			if (vm.count("profile"))  */
/* 				profile_flag = true; */
/* 			if (vm.count("alignment"))  */
/* 				align_flag = true; */
/* 			if (vm.count("verbose"))  */
/* 				verbose = true; */

/* 			if ( path_flag ) { */
/* 				if ( ! vm.count("input-list") && !vm.count("input-file") ) { */
/* 					std::cerr << "\nsequence file(s) must be given\n"; */
/* 					std::cerr << "\n" << USAGE_ASSEM << "\n\n"; */
/* 					std::cerr << visible; */
/* 					exit(1); */
/* 				} */
/* 			} */

/* 			if ( vm.count("input-list")) { */
/* 				file_list = vm["input-list"].as<std::string>(); */
/* 				loadFileList(file_list.c_str(), input_files); */
/* 			} */
/* 			if (vm.count("input-file")) { */
/* 				input_files = vm["input-file"].as< std::vector<std::string> >(); */
/* 			} */

/* 			checkKmerSize(k); */
/* 		} */
/* 		catch(std::exception& e) { */
/* 			std::cout << e.what() << "\n"; */
/* 			exit(1); */
/* 		}     */
/* 	} */

	void parse_cmd_args_assembly_long( int ac,
									   char **av,
									   Param &param)
	{
		try {
			std::string file_list;

			po::options_description generic("Generic options");
			generic.add_options()
				("version,v", "Print current version.")
				("help,h", "Produce help message.")
				;

			po::options_description graph("Graph options");
			graph.add_options()
				("kmer-size,k", po::value<int>(&param.kmer_size )->default_value(KMER_SIZE), "Size of k-mer (node)")
				("min-coverage,c", po::value<int>(&param.min_depth)->default_value(MIN_DEPTH), "Minimum kmer coverage for graph trimming")
				("no-trim", "Skip trimming graph")
				;

			po::options_description path("Path search options");
			path.add_options()
				("percentile,p", po::value<double>(&param.kmer_percentile)->default_value(KMER_PERCENTILE), "Seed k-mer: top x% percentile of kmers")
				("min-seed,s", po::value<int>(&param.min_seed)->default_value(SEED_START), "Minimum coverage of seed k-mer")
				("distance,d", po::value<int>(&param.back_trace)->default_value(BACK_TRACE), "Distance in path to find maximum read sharing node")
				("min-share,r", po::value<int>(&param.min_share)->default_value(SHARED_READS), "Minimum shared reads with neighboring node in path search")
				("min-length,n", po::value<int>(&param.min_length)->default_value(MIN_LENGTH), "Minimum path length")
				("no-path", "Skip path search")
				;

			po::options_description merge("Path merging/extension options");
			merge.add_options()
				("long-overlap,L", po::value<int>(&param.latch_length)->default_value(LATCH_LENGTH), "Minimum overlapping length to extend a sequence pair")
				("long-support,S", po::value<int>(&param.latch_support)->default_value(LATCH_SUPPORT), "MInimum number of supporting reads in ovelapping paths")
				("merge-score,m", po::value<double>(&param.merge_score)->default_value(MERGE_SCORE), "Score for merging sequences")
				("extend-score,e", po::value<double>(&param.latch_score)->default_value(LATCH_SCORE), "Score for extending sequences")
				("path-spur,u", po::value<int>(&param.path_spur)->default_value(PATH_MISMATCH_ALLOW), "Maximum length of spur (path to path)")
				("paired-score", po::value<double>(&param.paired_score)->default_value(PAIRED_SCORE), "Score for merging/extending paired paths")
				("filter-kmer", po::value<int>(&param.filter_kmer)->default_value(FILTER_KMER), "kmer size for sequence filtering")
				("no-merge", "Skip merging paths")
				;

			po::options_description read("Read recruitment options");
			read.add_options()
				("read-align,a", po::value<double>(&param.read_align_score)->default_value(READ_ALIGN_SCORE), "Minimum alignment score to be a part of a path")
				("read-ratio,y", po::value<double>(&param.read_align_ratio)->default_value(READ_ALIGN_RATIO), "Minimum length ratio of a read to be a part of a path")
				("read-spur,U", po::value<int>(&param.read_spur)->default_value(READ_MISMATCH_ALLOW), "Maximum length of spur (path to read)")
				("extend-reads", "Extend reads to a path")
				("no-recruit", "Skip read recruitment")
				;

			po::options_description pairend("Short path extension options");
			pairend.add_options()
				("extend-overlap,t", po::value<int>(&param.pairend_overlap)->default_value(PAIREND_OVERLAP), "Minimum overlapping length of short path extension")
				("extend-support,z", po::value<int>(&param.overlap_support)->default_value(OVERLAP_SUPPORT), "Minimum supporting paired end reads for short overlap extension")
				("extend-ratio,E", po::value<double>(&param.latch_ratio)->default_value(LATCH_RATIO), "Minimum percentage of read length to be extended to a path")
				("no-extend", "Skip path extension")
				("scaffold", "Scaffold paths")
				("correction", "Paired end read correction")
				//("tune", "Tune paths")
				;

			po::options_description other("Other program options");
			other.add_options()
				("med-coverage,C", po::value<int>(&param.med_depth)->default_value(MED_DEPTH), "Minimum median read coverage of amino acid in consensus")
				("base-depth,b", po::value<int>(&param.base_depth)->default_value(BASE_DEPTH), "Minimum base coverage of amino acid in a path")
				("filter-score,f", po::value<double>(&param.filter_score)->default_value(FILTER_SCORE), "Score for kmer filtering")
				("gap-open,g", po::value<int>(&param.gap_open)->default_value(GAP_OPEN), "Gap open penalty")
				("gap-ext,x", po::value<int>(&param.gap_ext)->default_value(GAP_EXT), "Gap extension penalty")
				("seed,R", po::value<unsigned>(&param.seed)->default_value(RAND_SEED), "random seed")
				("platform,T", po::value<int>(&param.platform)->default_value(PLATFORM), "Sequencing technology")
				("insert-size,w", po::value<int>(&param.insert_size)->default_value(INSERT_SIZE), "Inner insert size")
				("insert-sd,W", po::value<double>(&param.insert_sd)->default_value(INSERT_SD), "Insert std. deviation")
				("verbose,V", "verbose switch on")
				;

			po::options_description list("File options");
			list.add_options()
				("input-file,i", po::value< std::vector<std::string> >(), "A sequence file")
				("input-list,l", po::value<std::string>(&file_list), "A list of sequence files")
				("graph-file,G", po::value<std::string>(&param.graph_input), "Graph file")
				("index-file,I", po::value<std::string>(&param.index_file), "Inverted index file")
				("output-dir,o", po::value<std::string>(&param.out_dir)->default_value(OUTDIR), "Output directory")
				("file-suffix,F", po::value<std::string>(&param.file_suffix), "File suffix for loading pre-written dump")
				("dump-suffix,X", po::value<std::string>(&param.dump_suffix), "File suffix for writing dump")
				("debug-file,B", po::value<std::string>(&param.debug_file), "K-mer list for debugging purpose")
				("dump-flag,D", "Dump objects")
				("pair-end,M", "Paired end reads")
				("profile,P", "Generate a profile file")
				("alignment,A", "Generate an alignment file")
				("no-output", "Do not generate output files")
				;

			// Hidden options, will be allowed both on command line and
			// in config file, but will not be shown to the user.
			po::options_description hidden("Hidden options");
			hidden.add_options()
				("input-file(s)", po::value< std::vector<std::string> >(), "input file")
				;


			po::options_description cmdline_options;
			//cmdline_options.add(generic).add(parameter).add(list).add(hidden);
			cmdline_options.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list).add(hidden);

			po::options_description visible("Options");
			//visible.add(generic).add(parameter).add(list);
			visible.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list);

			po::positional_options_description p;
			p.add("input-file", -1);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).positional(p).run(), vm);
			notify(vm);


			if (vm.count("version")) {
				std::cout << VERSION << "\n";
				exit (1);
			}
			if (vm.count("help")) {
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("input-list") && !vm.count("input-file") ) {
				std::cerr << "\n[Error] sequence file(s) must be given\n";
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible;
				exit(1);
			}

			if ( !vm.count("file-suffix") && !vm.count("graph-file") ) {
				std::cerr << "\n[Error] K-mer to read mapping file must be provided by either directly or dump\n";
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible;
				exit(1);
			}
			if ( !vm.count("file-suffix") && !vm.count("index-file") ) {
				std::cerr << "\n[Error] K-mer to read mapping file must be provided by either directly or dump\n";
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible;
				exit(1);
			}

			if (vm.count("pair-end"))
				param.pair_flag = true;
			if (vm.count("no-trim"))
				param.trim_flag = false;
			if (vm.count("no-path"))
				param.path_flag = false;
			if (vm.count("no-merge"))
				param.merge_flag = false;
			if (vm.count("no-recruit"))
				param.recruit_flag = false;
			if (vm.count("no-extend"))
				param.extend_flag = false;
			if (vm.count("scaffold"))
				param.scaffold_flag = true;
			if (vm.count("correction"))
				param.correction_flag = true;
/* 			if (vm.count("tune")) */
/* 				param.tune_flag = true; */
			if (vm.count("no-output"))
				param.output_flag = false;
			if (vm.count("extend-reads"))
				param.extend_read_flag = true;
			if (vm.count("dump-flag"))
				param.dump_flag = true;
			if (vm.count("profile"))
				param.profile_flag = true;
			if (vm.count("alignment"))
				param.align_flag = true;
			if (vm.count("verbose"))
				param.verbose = true;

/* 			if ( param.path_flag ) { */
/* 				if ( ! vm.count("input-list") && !vm.count("input-file") ) { */
/* 					std::cerr << "\nsequence file(s) must be given\n"; */
/* 					std::cout << "\n" << USAGE_ASSEM << "\n\n"; */
/* 					std::cout << visible; */
/* 					exit(1); */
/* 				} */
/* 			} */

			if ( vm.count("input-list")) {
				file_list = vm["input-list"].as<std::string>();
				loadFileList(file_list.c_str(), param.input_files);
			}
			if (vm.count("input-file")) {
				param.input_files = vm["input-file"].as< std::vector<std::string> >();
			}
			
			checkKmerSize(param.kmer_size);
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}
		
	}
	
	void parse_cmd_args_assembly_short( int ac,
										char **av,
									   Param &param)
	{
		try {
			std::string file_list;

			po::options_description generic("");
			generic.add_options()
				//("version,v", "Current version.")
				("help,h", "This help message")

				("pair-end,M", "Paired-end reads [Input reads needs to be Interlaced]")
				("input-file,i", po::value< std::vector<std::string> >(), "Sequence file(s)")
				("input-list,l", po::value<std::string>(&file_list), "File of multiple sequence locations")
				("graph-file,G", po::value<std::string>(&param.graph_input), "Graph file")
				("index-file,I", po::value<std::string>(&param.index_file), "Inverted index file")
				("output-dir,o", po::value<std::string>(&param.out_dir)->default_value(OUTDIR), "Output directory")

				("kmer-size,k", po::value<int>(&param.kmer_size )->default_value(KMER_SIZE), "Size of k-mer (node)")
				("min-coverage,c", po::value<int>(&param.min_depth)->default_value(MIN_DEPTH), "Minimum kmer coverage for trimming graph")
				("percentile,p", po::value<double>(&param.kmer_percentile)->default_value(KMER_PERCENTILE), "Top x% percentile of seed kmers")
				("distance,d", po::value<int>(&param.back_trace)->default_value(BACK_TRACE), "Maximum distant node in a path to find same reads in a current node")
				//("min-share,r", po::value<int>(&param.min_share)->default_value(SHARED_READS), "Minimum number of same reads in a neighboring node\n")
				("min-share,r", po::value<int>(&param.min_share)->default_value(SHARED_READS), "Minimum node depth")
				("med-coverage,C", po::value<int>(&param.med_depth)->default_value(MED_DEPTH), "Median base depth in consensus")
				//("base-depth,b", po::value<int>(&param.base_depth)->default_value(BASE_DEPTH), "Minimum base coverage of amino acid in a path")

				("file-suffix,F", po::value<std::string>(&param.file_suffix), "File suffix for loading pre-written binrary objects")
				("dump-suffix,X", po::value<std::string>(&param.dump_suffix), "File suffix for writing binary objects")
				("dump-flag,D", "Dump objects")
				("profile,P", "Write a profile file")
				("alignment,A", "Write an alignment file")
				("verbose,V", "verbose switch on\n")
				("no-trim", "Skip trimming graph")
				("no-path", "Skip path search")
				("no-merge", "Skip merging paths")
				("no-recruit", "Skip read recruitment")
				("no-extend", "Skip path extension")
				//("no-correction", "Skip paired end read correction")

				//("gap-open,g", po::value<int>(&param.gap_open)->default_value(GAP_OPEN), "Gap open penalty")
				//("gap-ext,x", po::value<int>(&param.gap_ext)->default_value(GAP_EXT), "Gap extension penalty")

				//("no-output", "Do not generate output files")
				;

			// Hidden options, will be allowed both on command line and
			// in config file, but will not be shown to the user.
			po::options_description hidden("Hidden options");
			hidden.add_options()
				("input-file(s)", po::value< std::vector<std::string> >(), "input file")
				;


			po::options_description cmdline_options;
			//cmdline_options.add(generic).add(parameter).add(list).add(hidden);
			//cmdline_options.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list).add(hidden);
			cmdline_options.add(generic).add(hidden);

			po::options_description visible("Options");
			//visible.add(generic).add(parameter).add(list);
			//visible.add(generic).add(graph).add(path).add(merge).add(read).add(pairend).add(other).add(list);
			visible.add(generic);

			po::positional_options_description p;
			p.add("input-file", -1);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).positional(p).run(), vm);
			notify(vm);


			if (vm.count("version")) {
				std::cout << VERSION << "\n";
				exit (1);
			}
			if (vm.count("help")) {
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			//if ( param.path_flag ) {
				if ( ! vm.count("input-list") && !vm.count("input-file") ) {
					std::cerr << "\n[Error] sequence file(s) must be given\n";
					std::cout << "\n" << USAGE_ASSEM << "\n\n";
					std::cout << visible;
					exit(1);
				}
				//}

			if ( !vm.count("file-suffix") && !vm.count("graph-file") ) {
				std::cerr << "\n[Error] K-mer to read mapping file must be provided by either directly or dump\n";
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible;
				exit(1);
			}
			if ( !vm.count("file-suffix") && !vm.count("index-file") ) {
				std::cerr << "\n[Error] k-mer to read mapping file must be provided by either directly or dump\n";
				std::cout << "\n" << USAGE_ASSEM << "\n\n";
				std::cout << visible;
				exit(1);
			}

			if (vm.count("pair-end"))
				param.pair_flag = true;
			if (vm.count("no-trim"))
				param.trim_flag = false;
			if (vm.count("no-path"))
				param.path_flag = false;
			if (vm.count("no-merge"))
				param.merge_flag = false;
			if (vm.count("no-recruit"))
				param.recruit_flag = false;
			if (vm.count("no-extend"))
				param.extend_flag = false;
			if (vm.count("scaffold"))
				param.scaffold_flag = true;
/* 			if (vm.count("correction")) */
/* 				param.correction_flag = true; */
/* 			if (vm.count("no-correction")) */
/* 				param.correction_flag = false; */
/* 			if (vm.count("tune")) */
/* 				param.tune_flag = true; */
			if (vm.count("no-output"))
				param.output_flag = false;
			if (vm.count("extend-reads"))
				param.extend_read_flag = true;
			if (vm.count("dump-flag"))
				param.dump_flag = true;
			if (vm.count("profile"))
				param.profile_flag = true;
			if (vm.count("alignment"))
				param.align_flag = true;
			if (vm.count("verbose"))
				param.verbose = true;


			if ( vm.count("input-list")) {
				file_list = vm["input-list"].as<std::string>();
				loadFileList(file_list.c_str(), param.input_files);
			}
			if (vm.count("input-file")) {
				param.input_files = vm["input-file"].as< std::vector<std::string> >();
			}

			checkKmerSize(param.kmer_size);
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}

	}

	/**
	 * Parse command line arguements for assembler.
	 */
	void parse_cmd_args_assembly( int ac,
							  	  char **av,
							  	  Param &param)
	{

#if LONGKMER == 1
		parse_cmd_args_assembly_long(ac, av, param);
#else
		parse_cmd_args_assembly_short(ac, av, param);
#endif

	}

	/**
	 * Parse command line arguements for post processing.
	 */ 
	void parse_cmd_args_postspa( int ac, 
								 char **av, 
								 std::string &faa_file,
								 std::string &ffn_file,
								 std::string &place_file,
								 std::string &out_file )
	{
		try {

			po::options_description generic("Generic options");
			generic.add_options()
				("version,v", "print version string")
				("help,h", "produce help message")
				;

			po::options_description list("File options");
			list.add_options()
				("faa-file,a", po::value<std::string>(&faa_file), "ffa file")
				("ffn-file,d", po::value<std::string>(&ffn_file), "ffn file")
				("place-file,p", po::value<std::string>(&place_file), "place file")
				("out-file,o", po::value<std::string>(&out_file), "output file")
				;
		
			po::options_description cmdline_options;
			cmdline_options.add(generic).add(list);

			po::options_description visible("Options");
			visible.add(generic).add(list);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).run(), vm);
			notify(vm);

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_POST_SPA << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("faa-file") ||
				 ! vm.count("ffn-file") ||
				 ! vm.count("place-file") ||
				 ! vm.count("out-file")) {
				std::cerr << "[Error] .ffa .ffn .place output files must be given\n";
				std::cout << "\n" << USAGE_POST_SPA << "\n\n";
				std::cout << visible;
				exit(1);
			}
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    
	}


	/**
	 * Parse command line arguments for sequence pre-processing.
	 */ 
	void parse_cmd_metisin(int ac, 
						   char **av, 
						   int &k, 
						   std::string &graph_input,
						   std::string &metis_file,
						   std::string &nkmer_file )
	{
		try {
			/* Declare a group of options that will be allowed only on command line */
			po::options_description generic("Generic options");
			generic.add_options()
				("version,v", "print version string")
				("help,h", "produce help message")
				;
			
			po::options_description parameter("Program parameters");
			parameter.add_options()
				("kmer,k", po::value<int>(&k)->default_value(KMER_SIZE), "size of k-mer")
				;
			
			po::options_description list("File options");
			list.add_options()
				("graph-input,G", po::value<std::string>(&graph_input), "graph input name")
				("metis-file,m", po::value<std::string>(&metis_file), "metis input name")
				("nkmer-file,n", po::value<std::string>(&nkmer_file), "number to kmer mapping")
				;
		
			po::options_description cmdline_options;
			cmdline_options.add(generic).add(parameter).add(list);

			po::options_description visible("Options");
			visible.add(generic).add(parameter).add(list);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).run(), vm);
			notify(vm);

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_METISIN << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("graph-input") && 
				 ! vm.count("metis-file") && 
				 ! vm.count("nkmer-file")) {
				std::cerr << "[Error] File name(s) must be given\n";
				std::cout << "\n" << USAGE_METISIN << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
			/* Error check */
			checkKmerSize(k);
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    
	}


	/**
	 * Parse command line arguments for sequence pre-processing.
	 */ 
	void parse_cmd_devide(int ac, 
						  char **av, 
						  int &k, 
						  std::vector<std::string> &input_files, 
						  std::string &graph_input,
						  std::string &index_input,
						  std::string &metis_output,
						  unsigned &nparts,
						  std::string &kmermap_file,
						  std::string &output_dir )
	{
		try {
			std::string file_list;

			/* Declare a group of options that will be allowed only on command line */
			po::options_description generic("Generic options");
			generic.add_options()
				("version,v", "print version string")
				("help,h", "produce help message")
				;
			
			po::options_description parameter("Program parameters");
			parameter.add_options()
				("kmer,k", po::value<int>(&k)->default_value(KMER_SIZE), "size of k-mer")
				("part,p", po::value<unsigned>(&nparts), "size of partitions")
				;
			
			po::options_description list("File options");
			list.add_options()
				("input-list,l", po::value<std::string>(&file_list), "input file list")
				("input-file,i", po::value< std::vector<std::string> >(), "input file")
				("graph-input,G", po::value<std::string>(&graph_input), "graph input name")
				("index-input,I", po::value<std::string>(&index_input), "index input name")
				("metis-output,m", po::value<std::string>(&metis_output), "metis output name")
				("kmer-map,n", po::value<std::string>(&kmermap_file), "kmer to number mapping")
				("output-dir,o", po::value<std::string>(&output_dir)->default_value("."), "output directory")
				;

			/* Hidden options, will be allowed both on command line and
			 * in config file, but will not be shown to the user. */
			po::options_description hidden("Hidden options");
			hidden.add_options()
				("input-file(s)", po::value< std::vector<std::string> >(), "input file")
				;
		
			po::options_description cmdline_options;
			cmdline_options.add(generic).add(parameter).add(list).add(hidden);

			po::options_description visible("Options");
			visible.add(generic).add(parameter).add(list);

			po::positional_options_description p;
			p.add("input-file", -1);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).run(), vm);
			notify(vm);

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_DEVIDE << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("input-list") &&
				 ! vm.count("input-file") ) {
				std::cerr << "[Error] Input file(s) must be given\n";
				std::cout << "\n" << USAGE_DEVIDE << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
			if ( vm.count("input-list")) {
				file_list = vm["input-list"].as<std::string>();
				loadFileList(file_list.c_str(), input_files);
			}
			if (vm.count("input-file")) {
				input_files = vm["input-file"].as< std::vector<std::string> >();
			}
			
			if ( ! vm.count("graph-input") &&
				 ! vm.count("index-input") &&  
				 ! vm.count("metis-output") && 
				 ! vm.count("kmer-map")) {
				std::cerr << "[Error] File name(s) must be given\n";
				std::cout << "\n" << USAGE_DEVIDE << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
			/* Error check */
			checkKmerSize(k);
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    
	}

	/**
	 * Parse command line arguments for sequence pre-processing.
	 */ 
	void parse_cmd_merge(int ac, 
						 char **av, 
						 unsigned &nparts,
						 std::string &part_suffix,
						 std::string &output_dir )
	{
		try {

			/* Declare a group of options that will be allowed only on command line */
			po::options_description generic("Generic options");
			generic.add_options()
				("version,v", "print version string")
				("help,h", "produce help message")
				;
			
			po::options_description option("Program options");
			option.add_options()
				("part,p", po::value<unsigned>(&nparts), "size of partitions")
				("suffix,s", po::value<std::string>(&part_suffix)->default_value("path"), "suffix of path objects")
				("output-dir,o", po::value<std::string>(&output_dir)->default_value("."), "output directory")
				;

			po::options_description visible("Options");
			visible.add(generic).add(option);

			po::options_description cmdline_options;
			cmdline_options.add(generic).add(option);

			po::variables_map vm;
			store(po::command_line_parser(ac, av).
				  options(cmdline_options).run(), vm);
			notify(vm);

			if (vm.count("help")) {
				std::cout << "\n" << USAGE_DEVIDE << "\n\n";
				std::cout << visible << "\n";
				exit (1);
			}

			if (vm.count("version")) {
				std::cout << "Multiple sources example, version 1.0\n";
				exit (1);
			}

			if ( ! vm.count("part") ) {
				std::cout << "# partition is missing\n";
				std::cout << "\n" << USAGE_MERGE << "\n\n";
				std::cout << visible;
				exit(1);
			}
			
		}
		catch(std::exception& e) {
			std::cerr << e.what() << "\n";
			exit(1);
		}    
	}


}

#endif


