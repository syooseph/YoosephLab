/** 
 * \file      cmdargs_simple.h
 * \brief     Command line argument handler.
 * \details   This module parses command line arguments.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      Modified on Fri 2013-12-13 02:23:18 PM
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
//#include <boost/lexical_cast.hpp>
//#include <boost/format.hpp>

#include "core.h"

namespace po = boost::program_options;

/** 
 * \brief Command line argument handler.
 */
namespace args 
{
	//const std::string VERSION = "version 0.2 (with suffix array)";

	/**
	 * Usage message of post processing FragGeneScan program.
	 */
	const std::string USAGE_POSTFGS = 
		"Usage: postfgs [options]\n"
		" e.g.: postfgs -a fgs.faa -n fgs.ffn -o fgs.out reads1.fasta [reads2.fasta ...]\n";
	
	/**
	 * Usage message for no. of suffix array.
	 */
	const std::string USAGE_PART = 
		"Usage: part [options]\n"
		" e.g.: part -i input.faa\n";

	/**
	 * Usage message of preprocessing of sequences.
	 */
	const std::string USAGE_PREPARE = 
		"Usage: prespa [options]\n"
		" e.g.: prespa -k 6 -G graph -I index -i input.faa\n";
	/**
	 * Usage message of main assembler
	 */
	const std::string USAGE_ASSEM = 
		"Usage: spa [options] inputs\n"
		" e.g.: spa -k 6 -M -i input.faa\n";
	
	/**
	 * Usage message of post processing
	 */
	const std::string USAGE_POST_SPA = 
		"Usage: rtran [options] inputs\n"
		" e.g.: rtran -a fgs.faa -d fgs.ffn -p place.bin\n";


	/**
	 * Echo program command
	 */
	void printCommand(int ac, char **av, std::ostream &out);
	
	/**
	 * Check whether valid k-mer size is given.
	 */
	void checkKmerSize(int k);

	/**
	 * \param	listfile	a file with a list of sequence files.
	 * \param   input_files	names of input files
	 * \post	Vector of sequence file names.
	 */
	void loadFileList(const char *listfile, std::vector<std::string> &input_files) ;

	/**
	 * Parse command line arguments for FragGeneScan post processing
	 */ 
	void parse_cmd_args_postfgs(int ac, char **av, 
								std::vector<std::string> &input_files, 
								std::string &fgs_out,
								std::string &fgs_ffn,
								std::string &fgs_faa,
								std::string &out_dir,
								bool &verbose);

	/**
	 * Parse command line argument for suffix array partitions
	 */
	void parse_cmd_args_parts(int ac, 
							  char **av, 
							  std::vector<std::string> &input_files );
	
	/**
	 * Parse command line arguments for sequence pre-processing.
	 */ 
	void parse_cmd_args_prepare(int ac, 
								char **av, 
								int &k, 
								std::vector<std::string> &input_files, 
								std::string &outdir,
								int &npart,
								bool &reverse_sfa,
								bool &sfa_only,
								bool &verbose);

	/** 
	 * Parse and load SPA parameters
	 */
	void parse_cmd_args_assembly( int ac,
								  char **av,
								  Param &param);
	

	/**
	 * Parse command line arguements for post processing.
	 */ 
	void parse_cmd_args_postspa( int ac, 
								 char **av, 
								 std::string &faa_file,
								 std::string &ffn_file,
								 std::string &place_file,
								 std::string &out_file,
								 int &ncpus,
								 int &gap_ext,
								 int &gap_open,
								 //int &lower_band,
								 //int &upper_band,
								 double &band_ratio,
								 bool &banded,
								 bool &align,
								 bool &profile,
								 bool &verbose);

}

#endif


