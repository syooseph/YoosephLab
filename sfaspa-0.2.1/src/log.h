/** 
 * @file       log.h
 * @brief      logger for assembly stages
 * @detail     Detailed assembly status (counts and timing) for each status. 
 * @date       Modified on Thu 2014-02-27 05:12:52 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __SEARCH_LOG_H__
#define __SEARCH_LOG_H__

#include <stdio.h>
#include <stdlib.h>

/**
 * \brief Log for read validation during intial path serch
 */
struct ReadAlignLog
{
	double t_total;    ///< total read placement time
	double t_start;    ///< init
	double t_align;    ///< total alignment time
	double t_range;    ///< match reange search
	double t_place;    ///< placement
	double t_clean;    ///< trimming
	double t_count;    ///< count matrix of range from (k+2)-mer to nback-mer
	double t_reads;    ///< read alignment
	double t_update;   ///< update count 

	size_t ct_range;   ///< no. of read ranges to be placed
	size_t ct_place;   ///< no. of placed reads
	size_t ct_weak;    ///< no. of dropped reads because of weak score
	size_t ct_short;   ///< no. of dropped reads because of short length
	size_t ct_tiny;    ///< no. of dropped reads due to very short length
	
	size_t ct_align;   ///< no. of pairwise alignments between path and reads 
	size_t ct_substr;  ///< no. of exact match
	size_t ct_basecmp; ///< no. of base composition comparisons

	size_t ct_diff_aln; ///< path align region != read align region
	size_t ct_read_off; ///< read start off
	size_t ct_diff_len; ///< read align region != read length

	size_t ct_same_len;
	size_t ct_good_len;
	size_t ct_poor_len;

	/** Constructor */
	ReadAlignLog();
	
	/** Combine with other log */
	void add( ReadAlignLog &other );

	/** Initialize */
	void init();
};

/**
 * \brief Log for initial path search
 */
struct SearchLog
{
	size_t ct_reads, ct_assem;
	size_t ct_path_fail_save, ct_path_fail_srch, ct_path_fail_read;
	size_t ct_path_succ_save, ct_path_succ_srch, ct_path_succ_read;
	size_t ct_path_sim_node_succ, ct_path_mul_node_succ;
	size_t ct_path_sim_node_fail, ct_path_mul_node_fail;
	size_t ct_align_fail_save, ct_align_fail_srch, ct_align_fail_read;
	size_t ct_align_succ_save, ct_align_succ_srch, ct_align_succ_read;
	size_t ct_align_sim_node_succ, ct_align_mul_node_succ;
	size_t ct_align_sim_node_fail, ct_align_mul_node_fail;
	size_t ct_update_succ_save, ct_update_succ_srch, ct_update_succ_read;
	size_t ct_update_fail_save, ct_update_fail_srch, ct_update_fail_read;
	size_t ct_update_sim_node_succ, ct_update_mul_node_succ;
	size_t ct_update_sim_node_fail, ct_update_mul_node_fail;
	size_t ct_extend;
	size_t ct_graph;
	size_t ct_nback, ct_nback_hit, ct_nback_mis;
	size_t ct_iback, ct_iback_hit, ct_iback_mis;
	size_t ct_array;
	size_t ct_refine;
	size_t ct_endsrch;

	size_t srch_fail, align_fail, save_fail;
	size_t srch_fail_min, srch_fail_max, srch_fail_sum;
	size_t align_fail_min, align_fail_max, align_fail_sum;
	size_t save_fail_min, save_fail_max, save_fail_sum;
	size_t save_succ, save_succ_min, save_succ_max, save_succ_sum;

	//size_t ct_align_substr, ct_align_base, ct_align_align;

	size_t ct_check_fail_read, ct_check_fail_trim, ct_check_fail_len;

	//size_t ct_node_trace;


	size_t ct_path_cycle, ct_path_end, ct_path_weak;

	size_t ct_init;
	size_t ct_good_seed;
	size_t ct_good_path, ct_poor_path, ct_short_path;
	size_t ct_skip_seed; //, ct_skip_path, ct_skip_align;
	size_t ct_align, ct_good_align;
	size_t ct_update_succ, ct_update_fail;

	size_t ct_used_start_succ, ct_used_start_fail;
	size_t ct_used_end_succ, ct_used_end_fail;

	size_t ct_remove_edge, ct_remove_node;

	size_t ct_path_extend_left, ct_path_extend_right;


	double et_total;
	//double et_extend;
	double et_graph;
	double et_nback;
	double et_iback;
	double et_array;
	double et_refine;
	double et_endsrch;

	double et_skip, et_skip_seed;//, et_skip_path, et_skip_align;
	double et_good_path, et_poor_path, et_short_path;
	double et_init, et_init_seed, et_init_graph, et_init_prog;
	double et_seed;
	double et_path, et_path_init, et_path_path, et_path_post;
	double et_path_extend;
	double et_path_extend_left, et_path_extend_right;
	double et_path_pre, et_path_term, et_path_maxn, et_path_cycle, et_path_next;
	double et_node_ngb, et_node_max, et_node_cyc, et_node_trace;
	double et_node_trace_add;
	double et_node_array;
	double et_next_ngb, et_next_trace, et_next_max;
	double et_node_nback;

	double et_align;
	double et_update, et_update_path, et_update_depth, et_update_graph, et_update_nback, et_update_reads;//, et_update_pid;
	double et_update_seed;

	double et_lapse;

	double et_node_order;
	
	double et_used_start;
	double et_used_end;
	double et_update_rstarts;
	double et_update_rends;

	double et_remove_edge, et_remove_node;

	double et_used_start_srch, et_used_end_srch;
	
	double et_srch_fail, et_align_fail, et_save_fail, et_save_succ;
	SearchLog();
	void init();
	void printSummary();
	void printProgress(ReadAlignLog &read_log);
	void printOneLog();
};

/**
 * @brief Log for path clustering
 */
struct MergeLog
{
	double et_total;
	double et_index;
	//double et_kmers;
	double et_filter;
	double et_search;
	double et_path_index;
	double et_simpaths;
	double et_mtest;
	double et_section;
	double et_bases;
	double et_align;
	double et_update;

	size_t ct_lgap;
	size_t ct_egap;
	size_t ct_lgap_min, ct_lgap_max, ct_lgap_sum;
	size_t ct_egap_min, ct_egap_max, ct_egap_sum;
	
	size_t ct_allpaths;
	size_t ct_simpaths;
	size_t ct_bases;
	size_t ct_align;
	size_t ct_success;
	size_t ct_bases_succ;
	size_t ct_align_succ;
	size_t ct_align_succ_indel;

	size_t ct_align_succ_lfix;
	size_t ct_align_succ_efix;
	size_t ct_align_succ_lexp;
	size_t ct_align_succ_eexp;
	//size_t ct_align_succ_lgap;
	//size_t ct_align_succ_egap;
	size_t ct_align_fail_lgap;
	size_t ct_align_fail_egap;
	size_t ct_align_fail_short;
	size_t ct_align_fail_score;
	
	MergeLog();
	void printSummary();
};

/**
 * @brief Log for path extension
 */
struct LatchLog
{
	bool ljoin;
	bool pe_support;

	double et_total;
	double et_left;
	double et_right;
	double et_path;
	double et_path_filter, et_path_entry, et_path_sort, et_section; // long
	double et_reads; // short
	double et_latch;
	double et_latch_check;
	double et_bases;
	double et_align;
	double et_connect;

	size_t ct_path;
	size_t ct_bases;
	size_t ct_align;
	size_t ct_trial, ct_latch_succ;
	size_t ct_bases_succ, ct_align_succ;
	size_t ct_align_succ_indel;

	size_t ct_fail_pread;
	size_t ct_fail_short;

	LatchLog();
	void setLJoin(bool l);
	void setPESupport(bool p);
	void printSummary();
};

/**
 * @brief Log for path extension with bridging reads
 */
struct BridgeLog
{
	bool pe_support;

	double et_total;
	double et_latch;
	double et_latchable;

	double et_candidate;
	double et_reads;
	double et_drop;
	double et_evidence;
	double et_update;
	double et_sort;

	double et_read_list;
	double et_read_pivot;
	double et_read_match;
	double et_read_bad;
	double et_read_trim;

	double et_read_pair;

	size_t ct_trial;
	size_t ct_success;
	size_t ct_fail;
	size_t ct_candidates;
	size_t ct_cand_max, ct_cand_min;
	size_t ct_iter;
	
	size_t ct_read_all, ct_read_share, ct_read_align, ct_read_good;

	size_t ct_succ_first;
	size_t ct_fail_first;

	size_t ct_read_fail;
	size_t ct_drop_fail;
	size_t ct_evidence_fail;

    size_t ct_fail_pread;
	size_t ct_fail_short;

	BridgeLog();
	void setPESupport( bool p );
	void printSummary();
};


/**
 * @brief Log for recruitment of unassigned reads
 */
struct RecruitLog
{
	bool by_anchor;

	size_t ct_sub;
	size_t ct_aln;
	size_t ct_cmp;
	size_t ct_can;
	size_t sub_good;
	size_t aln_good;
	size_t ct_good;
	double et_sub;
	double et_pos;
	double et_aln;
	double et_anc;
	double et_can;

    double et_can_all, et_can_good;
    double et_sbjct_find, et_sbjct_get;

	RecruitLog();
	void print();
};

/** 
 * @brief Log for k-mer anchoring during read recruiment
 */
struct AnchorLog
{
	double et_init, et_init_read, et_init_path, et_init_insert;
	double et_find;
	double et_range;
	double et_map;
	double et_substr;
	double et_score;

	size_t ct_all;
	size_t ct_valid;
	size_t ct_exist;
	size_t ct_pass;

	AnchorLog();
	void print();
	void add(AnchorLog &other);
};


#endif
