/** 
 * @file       param.h
 * @brief      SPA program parameters
 * @date       Modified on Tue 2013-12-17 06:44:18 PM
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */


#ifndef SPAPARAM_H_
#define SPAPARAM_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include "default.h"

typedef std::vector<std::string> StrVec;

/**
 * \brief program parameters
 */
struct Param
{
	static int kmer_size;           ///< k-mer size
	static int gap_open;            ///< gap open penalty
	static int gap_ext;             ///< gap extension penalty
	//static int lower_band;
	//static int upper_band;
	static int min_depth;           ///< minimum k-mer coverage for graph trimming
	static int min_share;           ///< mininum shared read support in search
	static int min_seed;            ///< minimum k-mer coverage of path start
	static int min_length;          ///< minimum path length
	static int back_trace;          ///< maximum back trace string length of search same reads
	//static int dup_length;          ///< redundant path check length
	static int extend_length;       ///< minimum length of shared latch region
	//static int pairend_overlap;     ///< minimum overlap of paired reads latch
	static int overlap_support;     ///< minimum paired end reads for short overlap latching
	static int suffix_overlap;       ///< minimum short overlap with suffix array 
	static int bridge_overlap;      ///< length of bridging reads against an end of path
	static int merge_filter_kmer;   ///< clustering filter kmer
	static int merge_shared_nkmer;  ///< minimum number of shared kmer
	static int merge_expand_nbase;  ///< no. of bases to expand region to extract a substring of more room
	static int extend_anchor_kmer;  ///< anchor kmer size of path extension
	static int extend_anchor_mink;  ///< anchor minimum shared kmer
	//static int short_filter_kmer;   ///< short overlap filter kmer
	static int extend_filter_kmer;  ///< filter kmer size of path extension
	static int extend_off_nbase;    ///< no. of (indel) bases allowed to find extend region
	static int recruit_filter_kmer; ///< filter kmer size of read recruitment
	static int recruit_min_filter ;  ///< filter kmer size of read recruitment
	static int chunk_size;          ///< length of subpath to find read support using suffix array during path extension
	static int insert_size;         ///< sequence library inner insert size
	static int nparts;              ///< no. of suffix arrays
	static int ncpus;               ///< no. of cpus
	static int debug_id;            ///< path ID for debugging purpose
	static int line_length;         ///< MSA column length
	static int shift_length;        ///< Alignment wider boundary

	static int min_pair_reads;
	static int min_bridge_reads;

	static int  verbose ;              ///< verbose option

	static double read_align_score;     ///< alignment score for a read against a path
	static double read_align_ratio;     ///< alignment length for a read against a path
	static double insert_sd;            ///< standard deviation of insert size
	static double merge_filter_score;   ///< kmer filter score for clustering paths
	static double extend_filter_score;  ///< kmer filter score for path extension
	static double recruit_filter_score; ///< kmer filter score for read recruitment
	static double recruit_score;        ///< read recruitment alignment score
	static double recruit_ratio;        ///< read alignment length ratio
	static double merge_score;          ///< score cutoff for clustering paths
	static double extend_score;         ///< score cutoff for path extension
	static double kmer_percentile;      ///< percentage of total kmers being used as seeds
	static double band_ratio;           ///< percentage of length for banded length
	static double merge_short_ratio;    ///< shorter sequence length alignment ratio

	/* static std::string graph_input;     ///< k-mer list input file */
	/* static std::string index_file;      ///< k-mer to read mapping */
	static std::string out_dir;        ///< output directory
	static StrVec input_files;         ///< multiple sequence files

	static bool pair_flag ;            ///< paired end
	static bool trim_flag ;            ///< graph trimming
	static bool path_flag ;            ///< path finding flag
	static bool check_flag ;           ///< check paths flag
	static bool merge_flag ;           ///< clustering paths
	//static bool recluster_flag ;           ///< reclustering paths
	static bool extend_flag ;          ///< extend path pairs
	static bool place_flag;            ///< path/read placement
	static bool recruit_flag ;         ///< read recuritment
	static bool report_flag;           ///< final output
	static bool short_overlap_only;    ///< skip long overlap extension
	static bool long_overlap_only;     ///< skip short overlap extension
	static bool read_overlap_only;     ///< skip short overlap extension
	static bool read_bridge_extend; 
	static bool identity_flag;         ///< use identity score insteady of similarity score during sequence comparison
	static bool recruit_by_align;      ///< recruit reads by alignment instead of quick anchor base comparison
	static bool dump_flag ;            ///< dump binary objects

	static bool profile_flag ;         ///< profile file flag
	static bool align_flag ;           ///< alignment file flag
	//static bool uniq_seed;            ///< true: do not use kmers discovered in a path for seed
	static bool seed_reuse;            ///< false: do not use kmers discovered in a path for seed
	static bool percent_flag;          ///< true: use percentage of seeds for path search
	static bool summary_flag;          ///< true: show brief summary during assembly
	static bool par_search;            ///< true: best neighbor search in parallel
	static bool skip_fail;             ///< false: skip nodes from failed path
	static bool try_later;           
	//static bool use_all_fail;
	static bool strict_length;         ///< false: strict legnth condition
	static bool short_trim_len;        ///< false: do not allow short sequence after trimming
	static bool path_node_check;       ///< check existence of zero coverage node
	static bool extend_first;          ///< path extension before clustering : false
	static bool check_stop;           ///< stop conflict check during extension : false
	//static bool drop_dup_path;        ///< stop conflict check during extension : false
	static bool check_cycle_entry;        ///< stop conflict check during extension : false
	static bool ignore_stop;          ///< ignore stop codon while extending
	static bool output_all;           ///< generate outputs of intemediate stages.
	static bool replace_stop;         ///< replace stop codon with other amino acid
	static bool clip_stop;            ///< clip sequence away after stop codon 

	static bool merge_fixed_kmer;    ///< fixed size kmer for filtering
	//static bool use_section;
	//static bool use_preset_filter;    ///< use preset kmer filter value
	static bool align_base_first;    ///< perform base comparison

	static bool banded_align;        ///< baded alignment
	static bool release_path;
	static bool ljoin_pe_support;           ///< paire end read support
	static bool sjoin_pe_support;           ///< paire end read support
	static bool rjoin_pe_support;
	static bool debug_flag;           

	static bool extend_cluster;      ///< allow leading/trailing gaps for clustering?
	//static bool super_verbose ;              ///< verbose flag	

	static bool right_ext_only;
	static bool exact_match_only;

	Param();
	void print(std::ostream &out);
};

#endif /* SPAPAR_H_ */
