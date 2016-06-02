#include "param.h"

Param::Param(){}

bool Param::pair_flag          = false; 
bool Param::trim_flag          = TRIM_GRAPH; 
bool Param::path_flag          = true;  
//bool Param::check_flag         = true;  
bool Param::merge_flag         = true;  
// bool Param::recluster_flag     = true;  
bool Param::recruit_flag       = true;  
bool Param::extend_flag        = true; 
bool Param::place_flag         = true;
bool Param::report_flag        = true;
bool Param::dump_flag          = DUMP_BINARIES;
bool Param::profile_flag       = false;
bool Param::align_flag         = false;
//bool Param::uniq_seed          = UNIQ_SEED;
bool Param::seed_reuse         = SEED_REUSE;
bool Param::identity_flag      = false;
bool Param::short_overlap_only = false;
bool Param::long_overlap_only  = false;
bool Param::read_overlap_only  = false;
bool Param::read_bridge_extend = false;
bool Param::recruit_by_align   = false;
bool Param::percent_flag       = false;
//bool Param::summary_flag       = false;
bool Param::par_search         = true;
bool Param::skip_fail          = SKIP_FAILED_SEED_FIRST;
bool Param::try_later          = RUN_FAILED_SEED_LATER;
//bool Param::use_all_fail       = true;
bool Param::strict_length      = false;
bool Param::short_trim_len     = false;
bool Param::path_node_check    = true;
bool Param::extend_first       = false;
bool Param::check_stop         = false;
//bool Param::drop_dup_path      = true;
bool Param::check_cycle_entry  = false;
bool Param::ignore_stop        = false;
bool Param::replace_stop       = REPLACE_STOP_CODON;
bool Param::clip_stop          = CLIP_AFTER_STOP_CODON;
bool Param::output_all         = OUTPUT_ALL;
//bool Param::use_section        = false;
bool Param::merge_fixed_kmer   = true;
//bool Param::use_preset_filter  = false;
bool Param::align_base_first   = false;
bool Param::banded_align       = true;
bool Param::release_path       = false;
bool Param::ljoin_pe_support   = LONG_OVERLAP_PATH_PAIREND_SUPPORT;
bool Param::sjoin_pe_support   = SHORT_OVERLAP_PATH_PAIREND_SUPPORT;
bool Param::rjoin_pe_support   = READ_BRIDGE_PATH_PAIREND_SUPPORT;
bool Param::debug_flag         = false;
bool Param::extend_cluster     = false;
//bool Param::super_verbose      = false;
bool Param::right_ext_only     = false;
bool Param::exact_match_only   = true;


int Param::kmer_size           = KMER_SIZE;
int Param::gap_open            = GAP_OPEN;
int Param::gap_ext             = GAP_EXT;
// int Param::lower_band          = LOWER_BAND;
// int Param::upper_band          = UPPER_BAND;
int Param::min_depth           = MIN_DEPTH;
int Param::min_share           = SHARED_READS;
int Param::min_seed            = MIN_SEED;
int Param::min_length          = MIN_LENGTH;
int Param::back_trace          = BACK_TRACE;
int Param::extend_length       = EXTEND_LENGTH;
//int Param::pairend_overlap     = PAIREND_OVERLAP;
int Param::suffix_overlap      = SUFFIX_OVERLAP;
int Param::bridge_overlap      = BRIDGE_OVERLAP;
int Param::merge_filter_kmer   = MERGE_FILTER_KMER;
int Param::merge_shared_nkmer  = MERGE_SHARED_NKMER;
int Param::merge_expand_nbase  = MERGE_EXPAND_NBASE;
//int Param::short_filter_kmer   = SHORT_FILTER_KMER;
int Param::extend_anchor_kmer  = EXTEND_ANCHOR_KMER;
int Param::extend_anchor_mink  = EXTEND_ANCHOR_MINK;
int Param::extend_filter_kmer  = EXTEND_FILTER_KMER;
int Param::extend_off_nbase    = EXTEND_OFF_NBASE;
int Param::recruit_filter_kmer = RECRUIT_FILTER_KMER;
int Param::recruit_min_filter  = RECRUIT_MIN_FILTER;
int Param::insert_size         = INSERT_SIZE;
int Param::chunk_size          = CHUNK_SIZE;
int Param::nparts              = NUM_PARTS;
int Param::ncpus               = NUM_CPUS;
int Param::debug_id            = -1;
int Param::line_length         = LINE_LENGTH;
int Param::shift_length        = SHIFT_LENGTH;
int Param::min_pair_reads      = MIN_PAIR_READS;
int Param::min_bridge_reads    = MIN_BRIDGE_READS;
int Param::verbose             = VERBOSITY;

double Param::insert_sd            = INSERT_SD;
double Param::read_align_score     = READ_ALIGN_SCORE;
double Param::read_align_ratio     = READ_ALIGN_RATIO;
double Param::merge_filter_score   = MERGE_FILTER_SCORE;
double Param::extend_filter_score  = EXTEND_FILTER_SCORE;
double Param::recruit_filter_score = RECRUIT_FILTER_SCORE;
double Param::recruit_score        = RECRUIT_SCORE;
double Param::recruit_ratio        = RECRUIT_RATIO;
double Param::merge_score          = MERGE_SCORE;
double Param::extend_score         = EXTEND_SCORE;
double Param::kmer_percentile      = KMER_PERCENTILE;
double Param::band_ratio           = BAND_RATIO;
double Param::merge_short_ratio    = MERGE_SHORT_RATIO;

std::string Param::out_dir  = OUTDIR;
StrVec Param::input_files;  


void Param::print(std::ostream &out)
{
    out << std::string(Param::line_length, '=') << "\n";
    out << "Program summary\n";
    out << "---------------\n";
    out << "Input file: ";
    for ( StrVec::iterator it = input_files.begin(); it != input_files.end(); ++it )
        out << *it << " ";
    out << "\n";
    out << "Paired end reads: ";
    pair_flag ? out << "yes\n" : out << "no\n";
    char rpath[1024];
    realpath(out_dir.c_str(), rpath);
    out << "Working directory: " << rpath << "\n";
    out << "Kmer size: " <<  kmer_size << "\n";
    out << "No. of suffix arrays: " << nparts << "\n";
    out << "No. of cpus: " << ncpus << "\n";


    out << "Alignment setting\n";
    out << "- Banded alignment: ";
    if ( banded_align ) {
        out << "yes\n" ;
        out << "- Band ratio (w.r.t path length): " << band_ratio << "\n";
    } else out << "no\n";
    out << "- Substitution matrix: BLOSUM62\n";
    out << "- gap open penalty: " << gap_open << "\n";
    out << "- gap extension penalty: " << gap_ext << "\n";

    out << "Path search: ";
    path_flag ? out << "yes\n" : out << "no\n";
    if ( path_flag ) {
        out << "- Graph trimming: ";
        trim_flag ? out << "yes\n" : out << "no\n";
        if ( trim_flag ) 
            out << "- Min. vertex/edge coverage for graph trimming: " << min_depth << "\n";
        if ( !percent_flag ) 
            out << "- Min. seed coverage: " << min_seed << "\n";
        else
            out << "- Percentage of seed kmers: " << kmer_percentile << "%\n";
        out << "- Skip seed k-mers discovered in another path: ";
        !seed_reuse ? out << "yes\n" : out << "no\n";
        //uniq_seed ? out << "yes\n" : out << "no\n";
        out << "- Back trace string length: " << back_trace << "\n";
        out << "- Min. reads support between two neighboring nodes: " << min_share << "\n";
        out << "- Min. path length at path extraction stage: " << min_length << "\n";
    }

    out << "Extending path: ";
    extend_flag ? out << "yes\n" : out << "no\n";
    out << "- Extension prior to clustering: ";
    extend_first ? out << "yes\n" : out << "no\n";
    if ( extend_flag ) {
        out << "- Perform read bridging path extension: ";
        read_bridge_extend ? out << "yes\n" : out << "no\n";
        if ( ! short_overlap_only ) {
            out << "- Kmer filter size for long overlap extension: " << extend_filter_kmer << "\n";
            out << "- Minimun length of long overlapping region: " << extend_length << "\n";
        }
        if ( ! long_overlap_only )
            out << "- Minimun length of short overlapping region: " << suffix_overlap << "\n";
        //     out << "- Kmer filter size for short overlap extension: " << short_filter_kmer << "\n";

        out << "- Score for overlap filtering: " << extend_filter_score << "\n";
        out << "- Score for overlap extension: " << extend_score << "\n";
        out << "- Pair-end read support requirement for long overlap extension: ";
        ljoin_pe_support ? out << "yes\n" : out << "no\n";
        out << "- Pair-end read support requirement for short overlap extension: ";
        sjoin_pe_support ? out << "yes\n" : out << "no\n";
        if ( read_bridge_extend ) {
            out << "- Pair-end read support requirement for read bridging extension: ";
            rjoin_pe_support ? out << "yes\n" : out << "no\n";
        }
    }

    out << "Clustering path: ";
    merge_flag ? out << "yes\n" : out << "no\n";
    if ( merge_flag ) {
        out << "- Minimum shared kmers: " << merge_shared_nkmer << "\n";
        out << "- Initial kmer fiter size for merging: " << merge_filter_kmer << "\n";
        out << "- Score for kmer filtering: " << merge_filter_score << "\n";
        out << "- Score for path merging: " << merge_score << "\n";
    }

    out << "Read placement: ";
    place_flag ? out << "yes\n" : out << "no\n";
    out << "Recruiting reads: ";
    recruit_flag ? out << "yes\n" : out << "no\n";
    if ( recruit_flag ) {
        out << "- Kmer filter size for read recruitment: " << recruit_filter_kmer << "\n";
        out << "- Score for recruitment filtering: " << recruit_filter_score << "\n";
        out << "- Score for read recruitment: " << read_align_score << "\n";
    }

    out << "Dump binary output: ";
    dump_flag ? out << "yes\n" : out << "no\n";
    
    out << std::string(Param::line_length, '=') << "\n";
}

