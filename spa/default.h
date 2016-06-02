/** 
 * \file default.h
 * Default parameters.
 * This header set default program parameters.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \copyright J. Craig Venter Institute.
 */


//==============================================================================
// Default program parameter setting
//==============================================================================

#ifndef __DEFAULT_H__
#define __DEFAULT_H__

const std::string VERSION = "v0.001";

//--------------
// Graph options
//--------------
const int KMER_SIZE    = 6;   ///< default k-mer size
const int MIN_DEPTH    = 5;   ///< default min. depth for graph trimming.

//------------------
// Search parameters
//------------------
const double KMER_PERCENTILE = 50;  ///< top x% seed kmers
const int    SEED_START      = 100; ///< min read support of seed kmers
const int    BACK_TRACE      = 20;  ///< distance of tracing the same reads (Illumina sequencing)
const int    SHARED_READS    = 5;   ///< minimum read support to traverse neighboring nodes
const int	 MIN_LENGTH      = 6;   ///< default minimum path length

//------------------------
// Merge/extend parameters
//------------------------
const int LATCH_LENGTH        = 10;   ///< minimum common sequence length for sequence latching
const int LATCH_SUPPORT       = 5;    ///< minimum read support for connecting paths
const double MERGE_SCORE      = 0.9;  ///< bubble, spur, frayed rope
const double LATCH_SCORE      = 0.9;  ///< latch score
const double PAIRED_SCORE     = 0.8;  ///< merging/extending paired path
const int PATH_MISMATCH_ALLOW = 5;    ///< spur & frayed rope (path search)
const int FILTER_KMER         = 5;    ///< 5-mer

//---------------------------
// Read recruiment parameters
//---------------------------
const int READ_MISMATCH_ALLOW = 5;   ///< spur & frayed rope (read recruitment)

//-------------------------------------------
// Short ovelapping path extension parameters
//-------------------------------------------
const int PAIREND_OVERLAP = 5;   ///< minimum common sequence length in paired end reads support
const int OVERLAP_SUPPORT = 5;   ///< minimum paired end reads for short overlap latching
const double LATCH_RATIO  = 0.5; ///<  minimum percentage of read length to be extended a path

//--------------
// Other options
//--------------
const int MED_DEPTH            = 10;   ///< median read support for consensus
const int BASE_DEPTH           = 3;    ///< minimum base coverage
const double FILTER_SCORE      = 0.85; ///< quick k-mer filtering score
const double READ_ALIGN_RATIO  = 0.9;  ///< a portion of read aligned to a path
const double READ_ALIGN_SCORE  = 0.9;  ///< a score of read aligned to a path
const int GAP_OPEN             = -11;  ///< gap open penalty
const int GAP_EXT              = -1;   ///< gap extension penalty
const int PLATFORM             = 0;    ///< Illumina sequencing.
const int INSERT_SIZE          = 300;  ///< Typical insert size for Illumina platform.
const int INSERT_SD            = 50;   ///< Typical indert standard deviation of Illumina platform.
const unsigned RAND_SEED       = 1;    ///< default random seed
const std::string OUTDIR       = ".";  ///< output directory

#endif
