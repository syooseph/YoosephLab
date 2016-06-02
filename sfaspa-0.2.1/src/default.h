/** 
 * \file      default.h
 * \brief     Default parameters.
 * \details   This header set default program parameters.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      Modified on Wed 2013-10-23 05:50:39 PM
 * \copyright J. Craig Venter Institute.
 */


//==============================================================================
// Default program parameter setting
//==============================================================================

#ifndef __DEFAULT_H__
#define __DEFAULT_H__

const std::string VERSION = "SFA-SPA 0.2.1 - Suffix array based SPA (2015-02-10)";

//--------------
// Graph options
//--------------
const int KMER_SIZE    = 6;   ///< default k-mer size
const int MIN_DEPTH    = 2;   ///< default min. depth for graph trimming.
const bool TRIM_GRAPH  = false;

//------------------
// Search parameters
//------------------
const double KMER_PERCENTILE = 50;  ///< top x% seed kmers
const int    MIN_SEED        = 5;   ///< minimum seed depth
//const bool   UNIQ_SEED       = true; ///< do not use kmer seeds found in another path
const bool   SEED_REUSE       = false; ///< do not use kmer seeds found in another path
const int    BACK_TRACE      = 15;  ///< distance of tracing the same reads (Illumina sequencing)
const int    SHARED_READS    = 5;   ///< minimum read support to traverse neighboring nodes
const int	 MIN_LENGTH      = 8;   ///< default minimum path length
const bool   SKIP_FAILED_SEED_FIRST = true; ///< skip k-mer seeds from failed paths
const bool   RUN_FAILED_SEED_LATER = true; ///< use failed seeds in 2nd round
const bool   IGNORE_STOP_CODON = false;  ///< keep extract even if stop codon is encountered

//------------------------
// Clustering parameters
//------------------------
const double MERGE_SCORE         = 0.9;  ///< minimum path clustering alignment score
const int MERGE_FILTER_KMER      = 5;    ///< path clustering filter kmer
const int MERGE_SHARED_NKMER     = 2;    ///< minimum shared kmers
const int MERGE_EXPAND_NBASE     = 5;    ///< no. of addition bases to cover
const double MERGE_FILTER_SCORE  = 0.85; ///< path clustering kmer filter score
const double MERGE_SHORT_RATIO   = 0.9;

//--------------------------
// Path extension parameters
//--------------------------
const int EXTEND_FILTER_KMER     = 2;    ///< path extension filter kmer
const int EXTEND_ANCHOR_KMER     = 6;    ///< path extension anchor kmer
const int EXTEND_ANCHOR_MINK     = 1;    ///< path extension minimum shared kmers
//const int SHORT_FILTER_KMER      = 6;    ///< short overlap path exension filter kmer
const int EXTEND_OFF_NBASE       = 1;    ///< allowance of (indel) during extension region detection
const double EXTEND_SCORE        = 0.9;  ///< minimum path extension alignment score
const double EXTEND_FILTER_SCORE = 0.85; ///< path extension kmer filter score

const int EXTEND_LENGTH          = 10;   ///< minimum common sequence length for sequence latching
const int CHUNK_SIZE             = 20;   ///< sub-path length for suffix array search to find read support
//const int PAIREND_OVERLAP        = 5;    ///< minimum sequence overlap in short overlapping paths with pair-end read supports
const int SUFFIX_OVERLAP         = 4;    ///< minimum sequence overlap in short overlapping paths with read supports by suffix array
const int BRIDGE_OVERLAP         = 5;    ///< minimum length of read to end sequence overlap to find bridging reads

const int MIN_PAIR_READS         = 5;
const int MIN_BRIDGE_READS       = 5;

const bool LONG_OVERLAP_PATH_PAIREND_SUPPORT = false;
const bool SHORT_OVERLAP_PATH_PAIREND_SUPPORT = true;
const bool READ_BRIDGE_PATH_PAIREND_SUPPORT = true;

//---------------------------
// Read placement  parameters
//---------------------------
const double READ_ALIGN_RATIO  = 0.9;     ///< a portion of read aligned to a path
const double READ_ALIGN_SCORE  = 0.9;     ///< a score of read aligned to a path

//---------------------------
// Read recruitment  parameters
//---------------------------
const int RECRUIT_FILTER_KMER     = 5;    ///< read recruitment filter kmer
const int RECRUIT_MIN_FILTER      = 2;    ///< min. filter kmer 
const double RECRUIT_FILTER_SCORE = 0.85; ///< kmer filtering score for read recruitment
const double RECRUIT_RATIO        = 0.9;  ///< aligned portion of read to be recruited
const double RECRUIT_SCORE        = 0.9;  ///< alignment score of read to be recruited
const double BAND_RATIO           = 0.5;
const bool BAND_ALIGN             = false;

//--------------
// Other options
//--------------
const int GAP_OPEN             = -11;  ///< gap open penalty
const int GAP_EXT              = -1;   ///< gap extension penalty
/* const int LOWER_BAND           = -5; */
/* const int UPPER_BAND           = 5; */
const int PLATFORM             = 0;    ///< Illumina sequencing.
const int INSERT_SIZE          = 300;  ///< Typical insert size for Illumina platform.
const int INSERT_SD            = 50;   ///< Typical indert standard deviation of Illumina platform.
const std::string OUTDIR       = ".";  ///< output directory
const int NUM_PARTS            = 1;    ///< no. of suffix arrays
const int NUM_CPUS             = 1;    ///< no. of CPUs
const int LINE_LENGTH          = 70;   ///< MSA line length
const int SHIFT_LENGTH         = 5;    ///< MSA line length

const int VERBOSITY            = 0;

const bool OUTPUT_ALL          = false; ///< Generate outputs of all assembly stages
const bool DUMP_BINARIES       = false; ///< Dump binary outputs of each stage
const bool REPLACE_STOP_CODON  = true;  ///< When a stop codon is appeared in the middle of consensus, replace it with next best amino acid 
const bool CLIP_AFTER_STOP_CODON = true;///< Clip away sequences after stop codon for final SPA assembly output

#endif
