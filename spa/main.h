/** 
 * Main program.
 * This controls all the procedures of assembly. 
 * Followings are the main tasks.
 * -# Load reads
 * -# Build graph
 * -# Build inverted index
 * -# Path discovery
 * -# Merge/latch paths
 * -# Read recruitment
 * -# Path latch in paired end reads
 * -# Write result
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \pre       Generate graph input, inverted index beforehand.
 * \bug       None.
 * \warning   It may comsume a large amount of physical memory.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SPAMAIN_H__
#define __SPAMAIN_H__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <cmath>
#include "kmer.h"
#include "cmdargs.h"
#include "galignment.h"
#include "lalignment.h"
//#include "lgalignment.h"
//#include "assembly.h"
#include "graph.h"
#include "sequence.h"
#include "now.h"
#include "path.h"
#include "msa.h"
#include "eval.h"
#include "filter.h"
//#include "readpile.h"
//#include "pairend.h"
#include "coverage.h"
#include "container.h"
#include "param.h"
#include "assembly.h"
//#include "ungraph.h"


//--------------------
// function prototypes
//--------------------

void load( int argc, 
		   char *argv[], 
		   BitString *&bstrs, 
		   char *&strands, 
		   ReadId *&pairs, 
		   int &nreads, 
		   PathId *&used_reads, 
		   std::set<KmerType> &debug_kmers, 
		   Param &param );

void readParams(int argc, 
				char **argv, 
				Param &param);

void loadReads( BitString *&bstrs, 
				char *&strands, 
				ReadId *&pairs, 
				int &nreads, 
				Param &param );

//void build( DeBruijnGraph &graph, CoverageMap &kmer_coverage, VertexToKmerMap &vertex_map, Param &param );

void makeGraph( DeBruijnGraph &graph, 
                VertexToKmerMap &vertex_map,
                CoverageMap &kmer_coverage,
                PathToAlnMap &path2aln_map,
                PathId *used_reads,
                int nreads,
                Param &param);

void loadDump( DeBruijnGraph &graph, 
               VertexToKmerMap &vertex_map,
               CoverageMap &kmer_coverage,
               PathToAlnMap &path2aln_map,
               PathId *used_reads,
               int nreads,
               Param &param );




void loadDebugKmers( std::set<KmerType> &debug_kmers, Param &param );


void loadIndex(InvertedIndex &index, Param &param);
void initReadMembership( PathId *&used_reads, int nreads );
void initReadFlags( PathId *used_reads, int nreads );

int getUnusedReadCount( PathId *used_reads, int nreads );


void setStrands( char *strands, char **tags, int nreads );
void setBitstrs( BitString *bstrs, char **seqs, int nreads );
void setReadPairs( ReadId *pairs, 
                   char   **tags,
                   int    nreads );
void dropPairInfo(char **tags, int nreads);

void build( DeBruijnGraph &graph, 
			VertexToKmerMap &vertex_map,
			CoverageMap &kmer_coverage,
			PathToAlnMap &path2aln_map, 
			InvertedIndex &iindex, 
			int nreads, 
			PathId *used_reads, 
			Param &param );

void makeGraph( DeBruijnGraph &graph, 
                VertexToKmerMap &vertex_map,
                CoverageMap &kmer_coverage,
                PathToAlnMap &path2aln_map,
                PathId *used_reads,
                int nreads,
                Param &param);

void buildIndex(InvertedIndex &iindex, Param &param);

void assemble( DeBruijnGraph &graph, 
			   InvertedIndex &iindex, 
			   PathToAlnMap &path2aln_map, 
			   VertexToKmerMap &vertex_map,
			   CoverageMap &kmer_coverage,
			   BitString *bstrs, 
			   char *strands, 
			   ReadId *pairs, 
			   int nreads, 
			   PathId *used_reads,
			   std::set<KmerType> &debug_kmers,
			   Param &param );

void write( DeBruijnGraph &graph, 
			InvertedIndex &iindex,
			PathToAlnMap &path2aln_map, 
			VertexToKmerMap &vertex_map,
			CoverageMap &kmer_coverage,
			BitString *bstrs, 
			char *strands, 
			ReadId *pairs, 
			PathId *used_reads,
			int nreads, 
			Param &param );

void writeResult( PathToAlnMap &path2aln_map,
                  BitString *bstrs,
                  char *strands,
                  ReadId *pairs, 
                  PathId       *used_reads,
                  int           nreads,
                  Param &param );


//void writeConsensus( std::fstream &out, MSA &msa, int count );
void writeConsensus( std::fstream &out, SpaPath *spath, int count );

void writePlacement( std::fstream &out, SpaPath *spath, int count );

//void writeStatistic( std::fstream &out, MSA &msa, int count );
void writeStatistic( std::fstream &out, SpaPath *spath, int count );
//void writeAlignment( std::fstream &out, MSA &msa, SpaPath *spath, int count );
void writeAlignment( std::fstream &out, SpaPath *spath, int count, BitString *bstrs );
//void writeProfile( std::fstream &out, MSA &msa, int count );
void writeProfile( std::fstream &out, SpaPath *spath, int count );

void writeUnusedKmers( const char *file,
                       DeBruijnGraph &graph,
                       VertexToKmerMap &vertex_map,
                       InvertedIndex &iindex,
                       Param &param );

void writeUnusedReads( PathId *read_flag,
					   int nreads,
					   Param &param
                       );

void dumpObjects(DeBruijnGraph &graph, 
               VertexToKmerMap &vertex_map,
               CoverageMap &kmer_coverage,
               InvertedIndex &iindex, 
               PathToAlnMap &path2aln_map,
               Param &param );

void release( PathToAlnMap &path2aln_map,
              InvertedIndex &iindex,
              char *strands,
              BitString *bstrs,
              PathId *used_reads,
              ReadId *pairs,
              Param &param );

#endif

