/*==============================================================================
 * BWA reader
 * Author: Youngik Yang
 * Created: Thu 2015-07-30 09:15:34 PM
 * Last modified: Fri 2015-08-21 02:54:31 PM 
 *------------------------------------------------------------------------------
 * Input: BWA file (SAM)
 * Output: Matrix (row: reads, column: genomes)
 *==============================================================================*/

#ifndef __BWA_READER_H__
#define __BWA_READER_H__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <list>
#include <cassert>
#include <unistd.h>
#include <regex>
#include "pstream.h" // to read in samtools
#include "log.h"

// Avoid too small number for sbject size 
using sbjct_t = uint32_t; 
using query_t = uint32_t;
using ContigLengths  = std::unordered_map<std::string, size_t>; 
using SbjNumLengths  = std::unordered_map<sbjct_t, size_t>; 
using ContigToNumber = std::unordered_map<std::string, sbjct_t>;
using NumberToContig = std::unordered_map<sbjct_t, std::string>;
using SbjctList      = std::list<sbjct_t>;
using QueryMatches   = std::unordered_map<query_t, SbjctList>;
using SbjctCounts    = std::unordered_map<sbjct_t, size_t>; // match count
/* using QueryMatches   = std::map<query_t, SbjctList>; */
/* using SbjctCounts    = std::map<sbjct_t, size_t>; // match count */

// Alignment offset - BWA default edit distance
//const size_t ALLOW = 4;

class BwaReader
{
 private:
	std::string bwa_file;
	std::string outdir;
	bool verbose;
	bool pairend;
	bool set_contig;
	ContigLengths  conlens;
	SbjNumLengths  sbjlens;
	ContigToNumber con2num;
	NumberToContig num2con;
	QueryMatches   querys;
	SbjctCounts    sbjcts;
	size_t aln_offset;   // alignment offset for irregular alignment
	
 public:
	BwaReader();
	BwaReader(std::string &f, bool);
	void parse();
	void dump();
	void read();
	
	/* Settters */
	void setVerbose( bool &v )          { verbose = v; }
	void setDirectory( std::string &o ) { outdir = o; }

	/* Getter */
	SbjNumLengths*  getSbjctLengths() { return &sbjlens; }
	ContigToNumber* getContigNumbers() { return &con2num; }
	NumberToContig* getContigNames()   { return &num2con; }
	QueryMatches*   getHits()          { return &querys; }
	SbjctCounts*    getSbjHits()       { return &sbjcts; }
 private:
	void init();
	void readFile();
	void getHeader( std::string &line );
	void setContigs();	
	bool parseLine( std::string &line, query_t rnum, std::string &rname, size_t &edist );

	/** Check various type alignment 
	 * Clipped alignment (soft and hard)
	 * Spliced alignment
	 * Padded alignment
	 */
	bool goodCigar( std::string &cigar );
	
	bool properAlignment( std::string &cigar );

	/** Get the count of various type alignment */
	size_t getTypeCount( std::string &cigar, std::string &pattern );

	/** Extract edit distance */
	std::pair<bool,int> getEditDistance( std::string &col);

};

#endif
