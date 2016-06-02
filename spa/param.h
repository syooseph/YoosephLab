/*
 * spapar.h
 *
 *  Created on: Feb 10, 2012
 *      Author: padosori
 */

#ifndef SPAPARAM_H_
#define SPAPARAM_H_

#include <string>
#include <vector>

#include "default.h"

typedef std::vector<std::string> StrVec;

struct Param
{
	int kmer_size;           ///< k-mer size
	int gap_open;            ///< gap open penalty
	int gap_ext;             ///< gap extension penalty
	//int min_weight;          ///< minimum graph traverse weight
	int min_depth;           ///< minimum k-mer coverage for graph trimming
	int min_share;           ///< mininum shared read support in search
	int min_seed;            ///< minimum k-mer coverage of path start
	int med_depth;           ///< median read support for consensus
	int min_length;          ///< minimum path length
	int back_trace;          ///< maximum back trace of search same reads
	int base_depth;          ///< minimum base coverage
	int latch_length;        ///< minimum length of shared latch region
	int latch_support;       ///< minimum read support of latch region
	int pairend_overlap;     ///< minimum overlap of paired reads latch
	int overlap_support;     ///< minimum paired end reads for short overlap latching
	//int end_allow;           // maximum length of spur region
	int path_spur;           ///< maximum length of spur region during path search
	int read_spur;           ///< maximum length of spur region during read recruitment
	int filter_kmer;         ///< kmer size of quick filtering
	//int shift_offset;        // left/right shift offset of read placement
	int platform;            ///< sequencer
	int insert_size;         ///< inner insert size
	double read_align_score; ///< alignment score for a read against a path
	double read_align_ratio; ///< alignment length for a read against a path
	double insert_sd;        ///< standard deviation of insert size
	double filter_score;     ///< quick filtering score
	double merge_score;      ///< score cutoff for merging seeds
	double latch_score;      ///< score cutoff for merging seeds
	double latch_ratio;      ///< ratio of read length to be latched to a path
	double paired_score;
	//double latch_search;     ///< mininum rate of shared no. of kmers during latch paths
	//double recruit_search;   ///< mininum rate of shared no. of kmers during read recruitment
	double kmer_percentile;    ///< search start kmer's depth
	unsigned seed;           ///< seeding for random seed replacement
	std::string graph_input;   ///< k-mer list input file
	std::string index_file;  ///< k-mer to read mapping
	std::string out_dir;     ///< output directory
	std::string file_suffix; ///< dump file suffix for loading
	std::string dump_suffix; ///< dump file suffix for writing
	std::string debug_file;  ///< debug kmer files
	StrVec input_files;      ///< multiple sequence files
	bool pair_flag ;  ///< paired end
	bool trim_flag ;   ///< graph trimming
	bool path_flag ;   ///< path finding flag
	bool merge_flag ;  ///< latch paths
	bool recruit_flag ;///< read recuritment
	bool extend_flag ; ///< extend path pairs
	bool extend_read_flag ; ///< extend reads to a path end
	bool scaffold_flag ; ///< paired end path latch
	//bool tune_flag ; ///< paired end path latch
	bool correction_flag;
	bool output_flag; 
	//bool index_flag = true;  // debug kmers = index ?
	bool dump_flag ;  ///< dump graph
	bool verbose ;    ///< verbose flag
	//bool compressed ;  ///< gz file?
	bool profile_flag ; ///< profile file flag
	bool align_flag ;   ///< alignment file flag

	Param()
	{
		pair_flag        = false;  ///< paired end
		trim_flag        = true;   ///< graph trimming
		path_flag        = true;   ///< path finding flag
		merge_flag       = true;  ///< latch paths
		recruit_flag     = true;///< read recuritment
		extend_flag      = true; ///< paired end path latch
		scaffold_flag    = false; ///< scaffold paths
		extend_read_flag = true; ///< latch reads to path
		correction_flag  = true;
		output_flag      = true;
		dump_flag        = false;  ///< dump graph
		verbose          = false;    ///< verbose flag
		profile_flag     = false; ///< profile file flag
		align_flag       = false;   ///< alignment file flag
		
		kmer_size        = KMER_SIZE;
		gap_open         = GAP_OPEN;
		gap_ext          = GAP_EXT;
		min_depth        = MIN_DEPTH;
		min_share        = SHARED_READS;
		min_seed         = SEED_START;
		med_depth        = MED_DEPTH;
		min_length       = MIN_LENGTH;
		back_trace       = BACK_TRACE;
		base_depth       = BASE_DEPTH;
		latch_length     = LATCH_LENGTH;
		latch_support    = LATCH_SUPPORT;
		pairend_overlap  = PAIREND_OVERLAP;
		overlap_support  = OVERLAP_SUPPORT;
		path_spur        = PATH_MISMATCH_ALLOW;
		read_spur        = READ_MISMATCH_ALLOW;
		filter_kmer      = FILTER_KMER;
		
		platform         = PLATFORM;   
		insert_size      = INSERT_SIZE;
		insert_sd        = INSERT_SD;
		read_align_score = READ_ALIGN_SCORE;
		read_align_ratio = READ_ALIGN_RATIO;
		filter_score     = FILTER_SCORE;
		merge_score      = MERGE_SCORE;
		latch_score      = LATCH_SCORE;
		latch_ratio      = LATCH_RATIO;
		paired_score     = PAIRED_SCORE;
		//latch_search;     ///< mininum rate of shared no. of kmers during latch paths
		//recruit_search;   ///< mininum rate of shared no. of kmers during read recruitment
		kmer_percentile  = KMER_PERCENTILE;    ///< search start kmer's depth
		seed             = RAND_SEED;           ///< seeding for random seed replacement
		out_dir          = OUTDIR;
		
	}
};

#endif /* SPAPAR_H_ */
