//==============================================================================
// Sat 2011-05-28 12:41:27 AM
//==============================================================================

#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__

#include <iostream>
#include <fstream>
#include <unistd.h>
#include "kmer.h"
#include "cmdargs.h"
#include "timer.h"
#include "sequence.h"
//#include "math.h"
#include "InvertedIndex.h"
#include "coverage.h"


// locate invalid sequences
void trim( std::set<int> &bad_index, char **seqs, int size, int k);

// previous amino acids for each k-mers (graph input)
void setKmerLinks( LeftAAsCoverageMap &prev_AAs, CoverageMap &, char **seqs, int nreads, std::set<int> &, int k);

// remove graph input
void trash( LeftAAsCoverageMap &);

// unique kmers in a given sequence
KmerSet getKmerSet( char *cseq, int k );

// inverted index build
void buildIndex( InvertedIndex &InvInd, char **seqs, int nreads, CoverageMap &, std::set<int> &bad_index, int k );

// write graph input file
void writeGraphInput( const char* file, LeftAAsCoverageMap &prev_AAs, CoverageMap & );

void computeReadKmerStat(char **, InvertedIndex &, std::set<int> &, int , int , std::fstream &);

/* void coverageSummary( InvertedIndex &, std::ostream & ); */

#endif
