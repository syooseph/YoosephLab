/** 
 * @file reverser.h
 * @date Tue 2013-12-17 02:43:37 PM
 * @author Youngik Yang
 * @version 0.2
 * @brief Reverse translation
 * @details 
 * This generates DNA sequences from peptide sequences. 
 * @copyright J. Craig Venter Institute.
 */

#ifndef __REVERSER_H__
#define __REVERSER_H__

#include <sstream>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/copy.hpp>
#include "cmdargs_simple.h"
#include "codon.h"
#include "metapath.h"
#include "progress.h"

/** 
 * \brief Reverse translation from amino acids sequences to DNA sequences.
 */
class Reverser 
{
 private:
	std::string faa_file;     ///< short peptide reads file
	std::string ffn_file;     ///< short DNA reads file
	std::string aln_file;     ///< read placement file
	std::string out_dir;      ///< working directory

	int         ncpus;        ///< # cpus
	bool        align_flag;   ///< generate alignment?
	bool        profile_flag; ///< generate profile?
	bool        verbose;      ///< verbosity flag

	MetaPath    mpaths;
	char        **dnas;       ///< DNA reads
	char        **orfs;       ///< peptide reads
	int         norfs;        ///< count of sequence reads

    PathAligners NucAlns;     ///< an array of path alignments

	Codon        codon;       ///< codon table
	double       init_time;   ///< start time

	Progress     progress;    ///< progress

    std::fstream *outs;       ///< temporary output files
	std::fstream *logs;       ///< temporary log files
	std::fstream *alns;       ///< temporary alignment files
	std::fstream *pros;       ///< temporary profile files


	int gap_open;
	int gap_ext;

	bool banded_align;
	double band_ratio;
	//int lower_band;
	//int upper_band;

//==========================
//	Private member functions
//==========================
 private:
	/* /\** Runner - top level function *\/ */
	/* void run();  */

	/** Load read placements */
	void loadPlacements();

	/** Load DNA & peptide reads */
    void loadSequences();

	/**
	 * Trim DNA reads.
	 * FragGeneScan's DNA read may not be exact of tripples.
	 */
	void trimDnaReads();

	/** 
	 * Trim single DNA read 
	 */
	void trimSingleRead(int i);	

	/** Generate DNA sequences from peptide sequences */
	void generate();

	/** Combine multiple output files into one file */
	void combine();

	/** Generate DNA sequence from peptide sequence */
	void rtranslate( int thread,
					 int nid,
					 PathAligner &Aln,
					 PathAligner &NucAln );

	/** Open temporary files */
	void initFiles();

	/** Close temporary files */
	void closeFiles();
	
//=========================
//	Public member functions
//=========================
 public:

	/**
	 * Constructor
	 * @param ffn  a file name of nucleotide reads
	 * @param faa  a file name of peptide reads
	 * @param aln  a file name of SPA read placement 
	 * @param out  output directory
	 * @param cpus number of CPUs
	 * @param a    generate alignment?
	 * @param p    generate profile?
	 * @param v    verbosity
	 */
	Reverser( std::string ffn,
			  std::string faa,
			  std::string aln,
			  std::string out,
			  int cpus,
			  bool a,
			  bool p,
			  bool v );

	~Reverser();

	void run();

	void setGapPenalties(int gex, int gop) 
	{
		gap_ext = gex; gap_open = gop;
	}

	/* void setAlignmentBands(int l, int u) */
	/* { */
	/* 	lower_band = l; upper_band = u; */
	/* } */
	void setBandRatio( double b )
	{
		band_ratio = b;
	}

	void setBandedAlignment(bool b)
	{
		banded_align = b;
	}
};

#endif
