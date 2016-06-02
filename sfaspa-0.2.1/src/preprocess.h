/**
 * \file      preprocess.h
 * \brief     Preprocessing of SPA assembler
 * \details   This generates suffix array(s), graph input, and bad read indices if exists.
 * \author    Youngik Yang
 * \version   0.001
 * \date      Sat 2011-05-28 12:41:27 AM
 * \date      Modified on Tue 2013-12-17 06:31:29 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string>
//#include <cstring>
#include "kmer.h"
#include "cmdargs_simple.h"
#include "timer.h"
#include "sequence.h"
//#include "bitstr.h"
#include "path.h"
#include "gsa.h"

static bool VERBOSE = false;

void setVerbose( bool v ) { VERBOSE = v; }

// locate invalid sequences
void trim( std::set<int> &bad_index, char **seqs, int size, int k);

// minimun # of parts
int getMinNumParts ( char **seqs, int nreads );

void buildSuffixArray( char **seqs, int nreads, int nparts, std::string out_dir );
void buildReverseSuffixArray( char **seqs, int nreads, int nparts, std::string out_dir );

void reverseReads( char **rseqs,
				   char **seqs,
                   int nreads );

// previous amino acids for each k-mers (graph input)
void setKmerLinks( LeftAAsCoverageMap &prev_AAs, CoverageMap &, char **seqs, int nreads, std::set<int> &, int k);

// remove graph input
void trash( LeftAAsCoverageMap &);

// unique kmers in a given sequence
KmerSet getKmerSet( char *cseq, int k );

// write graph input file
void writeGraphInput( const char* file, LeftAAsCoverageMap &prev_AAs, CoverageMap & );

void writeBadIndex( const char* file, std::set<int> &bad_index );

#endif
