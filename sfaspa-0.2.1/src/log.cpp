#include "log.h"

ReadAlignLog::ReadAlignLog()
{
    init();
}

void ReadAlignLog::init()
{
    t_total = t_start = t_align = t_range = t_place = t_clean = t_count = t_reads = 0.0;		
    t_update = 0;

    ct_range = ct_place = ct_weak = ct_short = ct_tiny = 0;
    ct_align = ct_basecmp = ct_substr = 0;

    ct_diff_aln = ct_read_off = ct_diff_len = 0;

    ct_same_len = ct_good_len = ct_poor_len = 0;
}

void ReadAlignLog::add( ReadAlignLog &other ) 
{
    t_total += other.t_total;
    t_start += other.t_start;
    t_align += other.t_align;
    t_range += other.t_range;
    t_place += other.t_place;
    t_clean += other.t_clean;
    t_count += other.t_count;
    t_reads += other.t_reads;
    t_update += other.t_update;

    ct_range += other.ct_range;
    ct_place += other.ct_place;
    ct_weak  += other.ct_weak;		
    ct_short += other.ct_short;
    ct_tiny  += other.ct_tiny;
    
    ct_align += other.ct_align;
    ct_basecmp += other.ct_basecmp;
    ct_substr += other.ct_substr;

    ct_diff_aln += other.ct_diff_aln;
    ct_read_off += other.ct_read_off;
    ct_diff_len += other.ct_diff_len;

    ct_same_len += other.ct_same_len;
    ct_good_len += other.ct_good_len;
    ct_poor_len += other.ct_poor_len;

}


//////////////////////////////////////////////////////////////////////

SearchLog::SearchLog() 
{
	// ct_beg = 0; /* all search start */
	// ct_skp = 0; /* skipped start */
	// ct_val = 0; /* valid start */
	// ct_suc = 0; /* success */
	// ct_bad = 0; /* failure */
	// ct_sfa = 0; /* suffix array search */
	// ct_see = 0; /* # of extensions */
	// ct_gra = 0; /* gragh walk */
	// ct_sht = 0; /* short path */
	// //ct_dup = 0; /* redundant path */
	// //ct_sub = 0; /* sub path */
	// ct_wea = 0; /* weak support */
	// ct_ext = 0; /* extension fail */
	// ct_stp = 0; /* stop codon */
	// ct_end = 0; /* sink/source */
	// ct_nop = 0; /* bad start */
	// //ct_red = 0; /* total redundancy check */
	// //ct_rec = 0; /* total paths */
    // ct_align = 0;
    // ct_vap = 0;
    // ct_nba = 0;

    // ct_nrs = 0;
    // ct_ass = 0;

	// ct_ext_fail_srch = ct_ext_fail_read = ct_ext_fail_both = 0;
	// ct_chk_fail_srch = ct_chk_fail_read = ct_chk_fail_both = 0;
	// ct_ext_succ_srch = ct_ext_succ_read = ct_ext_succ_both = 0;
	// ct_chk_succ_srch = ct_chk_succ_read = ct_chk_succ_both = 0;

    // ct_ext_sim_node = ct_ext_mul_node = 0;
    // ct_chk_sim_node = ct_chk_mul_node = 0;

    // ct_skip = 0;

	// et_all = 0.0; /* total elapsed */
	// et_val = 0.0; /* valid start check time */
	// et_skp = 0.0; /* skip search time */
	// et_sfa = 0.0; /* suffix array search time */
	// et_gra = 0.0; /* graph traversal time */
	// et_ngb = 0.0; /* extract good neighbor */
	// //et_red = 0.0; /* duplicate path check time */
	// et_max = 0.0; /* next extenstion node determination time */
	// et_rec = 0.0; /* save path time */
	// et_cyc = 0.0; /* cycle check */
	// et_suc = 0.0; /* good path extraction */
	// et_bad = 0.0; /* bad path */
	// et_one = 0.0; /* one path search */
	// //et_str = 0.0;
	// //et_bak = 0.0;
	// et_end = 0.0;
	// et_stp = 0.0;
	// et_ngb = 0.0;
	// et_ext = 0.0;
	// et_nop = 0.0;
    // et_wea = 0.0;
    // et_add = 0.0;
	// et_ini = 0.0; // init before extensioin
	// et_pre = 0.0; // prepare 
	// et_pat = 0.0; // path time
	// et_pos = 0.0;
    // et_lft = 0.0;
    // et_rgt = 0.0;
    // et_ord = 0.0;
    // et_sum = 0.0;
    // et_upc = 0.0;
    // et_upg = 0.0;
    // et_upb = 0.0;
    // et_upr = 0.0;
    // et_chk = 0.0;
    // et_nba = 0.0;

    // et_upa = 0.0;
	// // lt_srh = 0.0; /* local serach time */
	// // lt_pre = 0.0; /* prepare */
	// // lt_pos = 0.0; /* path check */

    ct_reads = ct_assem = 0;

	ct_extend = 0;
	ct_graph = 0;
	ct_nback = ct_nback_hit = ct_nback_mis = 0;
    ct_iback = ct_iback_hit = ct_iback_mis = 0;
	ct_array = 0;
    ct_refine = 0;
    ct_endsrch = 0;
    
	srch_fail = align_fail = save_fail = 0;
	srch_fail_min = align_fail_min = save_fail_min = 1000000;
    srch_fail_max = align_fail_max = save_fail_max = 0;
	srch_fail_sum = align_fail_sum = save_fail_sum = 0;
	save_succ_min = 1000000;
    save_succ = save_succ_max = save_succ_sum = 0;
	ct_path_fail_save = ct_path_fail_srch = ct_path_fail_read = 0;
	ct_path_succ_save = ct_path_succ_srch = ct_path_succ_read = 0;
	ct_path_sim_node_succ = ct_path_mul_node_succ = 0;
	ct_path_sim_node_fail = ct_path_mul_node_fail = 0;
	ct_align_fail_save = ct_align_fail_srch = ct_align_fail_read = 0;
	ct_align_succ_save = ct_align_succ_srch = ct_align_succ_read = 0;
	ct_align_sim_node_succ = ct_align_mul_node_succ = 0;
	ct_align_sim_node_fail = ct_align_mul_node_fail = 0;
	ct_update_fail_save = ct_update_fail_srch = ct_update_fail_read = 0;
	ct_update_succ_save = ct_update_succ_srch = ct_update_succ_read = 0;
	ct_update_sim_node_succ = ct_update_mul_node_succ = 0;
	ct_update_sim_node_fail = ct_update_mul_node_fail = 0;

	//ct_align_substr = ct_align_base = ct_align_align = 0;

    ct_check_fail_read = ct_check_fail_trim = ct_check_fail_len = 0;

    //ct_node_trace = 0;

	ct_init = 0;
    ct_good_seed = 0;
	ct_good_path = ct_poor_path = ct_short_path = 0;
	//ct_skip = 
    ct_skip_seed = 0; //ct_skip_path = ct_skip_align = 0;
    ct_align = ct_good_align = 0;
    ct_update_succ = ct_update_fail = 0;

    ct_path_cycle = ct_path_end = ct_path_weak = 0;

    ct_used_start_succ = ct_used_start_fail = 0;
    ct_used_end_succ = ct_used_end_fail = 0;

    ct_remove_node = ct_remove_edge = 0;

    ct_path_extend_right = ct_path_extend_left = 0;

    et_total = 0.0;
	//et_extend = 0.0;
	et_graph = 0.0;
	et_nback = 0.0;
	et_iback = 0.0;
	et_array = 0.0;
	et_refine = 0.0;
	et_endsrch = 0.0;

	et_skip = et_skip_seed = 0.0; //= et_skip_path = et_skip_align = 0.0;
	et_good_path = et_poor_path = et_short_path = 0.0;
	et_init = 0.0;
    et_init_seed = et_init_graph = et_init_prog = 0.0;
	et_seed = 0.0;
    et_node_order = 0.0;
	et_path = et_path_init = et_path_path = et_path_post = 0.0;
    et_path_extend = 0.0;
    et_path_extend_left = et_path_extend_right = 0.0;
    et_node_ngb = et_node_max = et_node_cyc = et_node_trace = 0.0;
    et_node_trace_add = 0;
    et_path_pre = et_path_term = et_path_maxn = et_path_cycle = et_path_next = 0.0;

    et_next_ngb = et_next_trace = et_next_max = 0.0;
    et_node_array = 0.0;
    et_node_nback = 0.0;

	et_align = 0.0;
	et_update = et_update_path = et_update_depth = et_update_graph = et_update_nback = et_update_reads = et_update_seed = 0.0;
	et_lapse = 0.0;

    et_used_start = et_used_end = 0.0;
    et_update_rstarts = et_update_rends = 0.0;

    et_remove_node = et_remove_edge = 0;
    
    et_used_start_srch = et_used_end_srch = 0.0;
	et_srch_fail = et_align_fail =  et_save_fail = et_save_succ = 0.0;
}

void SearchLog::init()
{
    //lt_srh = lt_pre = lt_pos = 0.0;
}

void SearchLog::printSummary()
{
    printf( "# Attempt:%zu\n", ct_init );
    printf( "# Skipped:%zu\n", ct_skip_seed );
    printf( "# Extract:%zu\n", ct_good_seed);
    printf( "# Success:%zu \n", ct_good_path );
    printf( "# Failure:%zu\n", ct_poor_path );// (dup:%zu, short:%zu)\n", ct_bad, ct_sub, ct_sht);
    printf( "# Invalid paths:%zu\n", ct_good_path-ct_good_align );
    //printf( "# Outputs:%zu\n", ct_good_align);
    printf( "# Outputs:%zu\n", ct_update_succ);
    printf( "# Recruit:%zu/#Reads:%zu (%.2f%%)\n", ct_assem, ct_reads, 100*(double)ct_assem/ct_reads);
    // std::cout << "\nExtension stop criteria\n";
    // std::cout << "\t#Extension:" << ct_see << "\n";
    // std::cout << "\tRedundant path:" << ct_dup << "\n";
    // std::cout << "\tStop codon:"    << ct_stp << "\n";
    // std::cout << "\tSink/Source:"    << ct_end << "\n";
    // std::cout << "\tWeak support:"   << ct_ext << "\n";
    // std::cout << "\tBad start:"      << ct_nop << "\n";    
    printf( "\nElapsed time:%.2f sec\n", et_total );
    //printf( "N-back query:%.2f sec, # search:%zu\n", et_nback, ct_nback );
    printf( "Suffix array:%.2f sec, # search:%zu\n", et_array, ct_array );
    printf( "Refine array:%.2f sec, # search:%zu\n", et_refine, ct_refine );
    printf( "End supports:%.2f sec, # search:%zu\n", et_endsrch, ct_endsrch );
    printf( "Graph search:%.2f sec, # search:%zu\n", et_graph, ct_graph );
    //printf( "Redundany check:%.2f sec\t# search:%zu\n", et_red, ct_red );
}

void SearchLog::printProgress(ReadAlignLog &read_log)
{
    printf( "# Attempt:%zu\n", ct_init );                    
    printf( "# Skipped:%zu (%.2f sec)\n", ct_skip_seed, et_skip_seed );
    //printf( "# Extract:%zu (%.2f sec)\n", ct_good_seed, et_path );
    printf( "# Extract:%zu (%.2f sec: [init:%.2f, path:%.2f, post:%.2f])\n", ct_good_seed, et_path, et_path_init, et_path_path, et_path_post );
    //printf( "# Path extraction:%.2f sec (node-traveral:%.2f max-node:%.2f [trace-str:%.2f {init:%.2f add:%.2f str:%.2f} order:%.2f suffix:%.2f count:%.2f nback:%.2f] repeat-check:%.2f)\n", et_path_path, et_node_ngb, et_node_max, et_node_trace, et_node_nback, et_node_trace_add, (et_node_trace-et_node_trace_add-et_node_nback), et_node_order, et_node_array, et_iback, et_nback, et_path_cycle);//et_node_cyc );
    printf( "# Path extraction:%.2f sec (node-traveral:%.2f, max-search:%.2f, repeat-check:%.2f)\n", et_path_path, et_node_ngb, et_node_max, et_path_cycle);//et_node_cyc );
    printf( "# Extension:%.2f sec (init:%.2f, end-node:%.2f, neighbor:%.2f, repeat:%.2f, next:%.2f)\n", et_path_extend, et_path_pre, et_path_term, et_path_maxn, et_path_cycle, et_path_next );
    printf( "- Left-extension:%zu (%.2f sec), right-extension:%zu (%.2f sec)\n", ct_path_extend_left, et_path_extend_left, ct_path_extend_right, et_path_extend_right );
    printf( "- Max node search:%.2f sec (extract neighboring nodes:%.2f, max node:%.2f)\n", et_path_maxn, et_next_ngb, et_next_max );
    //printf( "# Trace string:%zu (%.2f)\n", ct_node_trace, et_node_trace);
    printf( "# Termination conditions: (path-ends:%zu, weak-node:%zu, cycle:%zu)\n", ct_path_end, ct_path_weak, ct_path_cycle);
    printf( "# Success:%zu (%.2f sec)\n", ct_good_path, et_good_path );
    printf( "- Seeds from failed paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", 
            ct_path_succ_srch, ct_path_succ_read, ct_path_succ_save );
    printf( "- Node type (simple-node:%zu\tcomplex-node:%zu)\n", ct_path_sim_node_succ, ct_path_mul_node_succ);
    printf( "# Failure:%zu (%.2f sec)\n", ct_poor_path, et_poor_path);// (dup:%zu, short:%zu)\n", ct_bad, ct_sub, ct_sht);
    if ( srch_fail ) 
        printf("- Max length:%zu, Min length:%zu, Avg length:%.2f\n", srch_fail_max, srch_fail_min, srch_fail_sum/(double)srch_fail);
    printf( "- Seeds from bad paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", 
            ct_path_fail_srch, ct_path_fail_read, ct_path_fail_save );
    printf( "- Node type (simple-node:%zu, complex-node:%zu)\n", ct_path_sim_node_fail, ct_path_mul_node_fail);
    //printf( "# Path check:%.2f sec (align-total:%.2f\tinit:%.2f align:%.2f [reads:%.2f clear:%.2f range:%.2f place:%.2f dedup:%.2f] trim:%.2f)\n", et_chk, read_log.t_total, read_log.t_start, read_log.t_align, read_log.t_reads, read_log.t_clear, read_log.t_range, read_log.t_place, read_log.t_dedup, read_log.t_clean );
    //printf( "# Path check:%.2f sec (range:%.2f place:%.2f dedup:%.2f trim:%.2f)\n", et_chk, read_log.t_range, read_log.t_place, read_log.t_dedup, read_log.t_clean );
    printf( "# Path validity check:%.2f sec (range:%.2f, place:%.2f, count-matrix:%.2f, trim:%.2f, bound-upate:%.2f )\n", et_align, read_log.t_range, read_log.t_place, read_log.t_count, read_log.t_clean, read_log.t_update );
    //printf( "# Node type (simple-node:%zu\tcomplex-node:%zu)\n", ct_align_sim_node_succ, ct_align_mul_node_succ);
    //printf( "# Placement type (substring:%zu)\n", read_log.ct_substr );
    printf( "# Valid paths:%zu\n", ct_good_align);
    printf( "# Placement type (substring:%zu, base-cmp:%zu, alignment:%zu)\n", read_log.ct_substr, read_log.ct_basecmp, read_log.ct_align );
    printf( "# Read type (full length:%zu, good length:%zu, poor length:%zu)\n",
            read_log.ct_same_len, read_log.ct_good_len, read_log.ct_poor_len );
    printf( "- Seeds from bad paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", ct_align_succ_srch, ct_align_succ_read, ct_align_succ_save );
    printf( "- Node type (simple-node:%zu, complex-node:%zu)\n", ct_align_sim_node_succ, ct_align_mul_node_succ);
    printf( "# Invalid paths:%zu\n", ct_good_path-ct_good_align );
    if ( align_fail ) 
        printf("- Max length:%zu, Min length:%zu, Avg length:%.2f\n", align_fail_max, align_fail_min, align_fail_sum/(double)align_fail);
    printf( "- Seeds from bad paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", ct_align_fail_srch, ct_align_fail_read, ct_align_fail_save );
    printf( "- Node type (simple-node:%zu, complex-node:%zu)\n", ct_align_sim_node_fail, ct_align_mul_node_fail);
    printf( "# Align-fail:%zu, Trim-fail:%zu, Length-fail:%zu\n", ct_check_fail_read, ct_check_fail_trim, ct_check_fail_len );
    //printf( "# Ranges:%zu Placements:%zu Diff-length:%zu Short-align-ratio:%zu Tiny-align-length:%zu Weak-similarity:%zu\n", read_log.ct_range, read_log.ct_place, read_log.ct_diff, read_log.ct_short, read_log.ct_tiny, read_log.ct_weak);
    //printf( "# Ranges:%zu Placements:%zu Short-align-ratio:%zu Tiny-align-length:%zu Weak-similarity:%zu\n", read_log.ct_range, read_log.ct_place, read_log.ct_short, read_log.ct_tiny, read_log.ct_weak);
    printf( "# Ranges:%zu, Placements:%zu, Diff-align-length:%zu, Read-start-off:%zu, Read-diff-len:%zu\n", read_log.ct_range, read_log.ct_place, read_log.ct_diff_aln, read_log.ct_read_off, read_log.ct_diff_len);
    printf( "# Saved paths:%zu\n", ct_update_succ );
    if ( save_succ ) 
        printf("- Max length:%zu, Min length:%zu, Sum length:%zu, Avg length:%.2f\n", save_succ_max, save_succ_min, save_succ_sum, save_succ_sum/(double)save_succ);
    printf( "- Seeds from bad paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", ct_update_succ_srch, ct_update_succ_read, ct_update_succ_save );
    printf( "- Node type (simple-node:%zu, complex-node:%zu)\n", ct_update_sim_node_succ, ct_update_mul_node_succ);
    printf( "# Failed paths:%zu\n", ct_update_fail );
    if ( save_fail ) 
        printf("- Max length:%zu, Min length:%zu, Avg length:%.2f\n", save_fail_max, save_fail_min, save_fail_sum/(double)save_fail);
    printf( "- Seeds from bad paths (extract-failed:%zu, align-failed:%zu, update-failed:%zu)\n", ct_update_fail_srch, ct_update_fail_read, ct_update_fail_save );
    printf( "- Node type (simple-node:%zu, complex-node:%zu)\n", ct_update_sim_node_fail, ct_update_mul_node_fail);
    //printf( "# Update pid:%.2f sec\n", et_update_pid );
    printf( "# Update path:%.2f sec\n", et_update_path );
    printf( "# Update coverage:%.2f sec\n", et_update_depth );
    printf( "# Update graph:%.2f sec\n", et_update_graph );
    //printf( "# Update nback:%.2f sec\n", et_update_nback );
    printf( "# Update read-starts:%.2f sec\n", et_update_rstarts);
    printf( "# Update read-ends:%.2f sec\n", et_update_rends);
    printf( "# Update reads:%.2f sec\n", et_update_reads );
    printf( "# Update seeds:%.2f sec\n", et_update_seed );
    printf( "# Trimming:%.2f sec (nodes:%zu [%.2f sec], edges:%zu [%.2f sec])\n", et_remove_node+et_remove_edge, ct_remove_node, et_remove_node, ct_remove_edge, et_remove_edge );
    printf( "# Path extraction (Extraction fail:%.2f sec, Placement fail:%.2f sec, Coverage fail:%.2f sec, Good path:%.2f sec)\n", et_srch_fail, et_align_fail, et_save_fail, et_save_succ );
    printf( "# Path search time summary (elapsed:%.2f init:%.2f seed-validity:%.2f path-extraction:%.2f path-check:%.2f update:%.2f skipped:%.2f loop-lapse:%.2f)\n", et_total, et_init, et_seed, et_path, et_align, et_update, et_skip, et_lapse );
    printf( "- Path search init details (elapsed:%.2f seed:%.2f graph:%.2f progress:%.2f)\n", et_init, et_init_seed, et_init_graph, et_init_prog);
    //printf( "# Reduandant:(%zu)\t", ct_dup );
    // printf( "# Sub path:%zu\t", ct_sub);
    // printf( "# Short Path:%zu\n", ct_sht );
    //printf( "# Extension:%zu\n", ct_extend );
    printf( "- Nodes:%zu (%.2f sec)\n", ct_graph, et_graph);
    printf( "- SFA query:%zu (%.2f sec)\n", ct_array, et_array);
    printf( "- Refine array:%zu (%.2f sec)\n", ct_refine, et_refine );
    printf( "- End supports:%zu (%.2f sec)\n", ct_endsrch, et_endsrch );
    //printf( "# N-back query:%zu (%.2f sec) [#hit:%zu #miss:%zu]\n", ct_nback, et_nback, ct_nback_hit, ct_nback_mis );
    printf( "- R-start query:%zu (%.2f sec) [#hit:%zu, #miss:%zu]\n", (ct_used_start_succ+ct_used_start_fail), et_used_start, ct_used_start_succ, ct_used_start_fail);
    printf( "- R-end query:%zu (%.2f sec) [#hit:%zu, #miss:%zu]\n", (ct_used_end_succ+ct_used_end_fail), et_used_end, ct_used_end_succ, ct_used_end_fail);
    printf( "- S-back query:%zu (%.2f sec) [#hit:%zu, #miss:%zu]\n", ct_iback, et_iback, ct_iback_hit, ct_iback_mis );
    printf( "# Recruit:%zu/#Reads:%zu (%.2f%%)\n", ct_assem, ct_reads, 100*(double)ct_assem/ct_reads);
    //printf( "# Redundany:%zu (%.2f sec)\n", ct_red, et_red );



    // printf( "Other timing detail\n");
    // printf( "Extraction preparation:%.2f\n", et_pre);
    // printf( "Path extraction:%.2f\n", et_pat);
    // printf( "\tLeft:%.2f\tRight:%.2f\n", et_lft, et_rgt);
    // printf( "\tInit:%.2f\n", et_ini);
    // printf( "\tNeighbor search:%.2f\n", et_ngb );
    // printf( "\tBest next node:%.2f\n", et_max );
    // printf( "\tSort:%.2f\tSum:%.2f\n", et_ord, et_sum );
    // printf( "\tTrace string:%.2f\n", et_str );
    // printf( "\tBack tracking:%.2f\n", et_bak );
    // printf( "\tPath end:%.2f\n", et_end );
    // printf( "\tStop codon:%.2f\n", et_stp );
    // printf( "\tWeak extension:%.2f\n", et_wea );
    // printf( "\tBad start:%.2f\n", et_nop );
    // printf( "\tCycle handling:%.2f\n", et_cyc );
    // printf( "\tOther:%.2f\n", et_add );
    // printf( "Path validation:%.2f\n", et_pos);
}

void SearchLog::printOneLog()
{
    printf( "#iterations:%zu\t", ct_init );
    printf( "search:%.4f sec\t", et_path );
    printf( "check:%.4f sec\t", et_align );
    printf( "update:%.4f\t", et_update );
    printf( "total:%.4f sec\t", et_total );
    printf( "average:%.4f sec\n", et_total/ct_init );
 //    printf( "success:%.4f sec\t", et_suc );
 //    printf( "failure:%.4f sec\t", et_bad );
    
 // skiped:%.2f\tneighbor:%.2f\tmax-node:%.2f\tredundant:%.2f\tsave:%.2f\tgood:%.2f\tbad:%.2f\n",
 //                    ct_beg, mytime()-tic, mytime()-lt0, (mytime()-lt0)/ct_beg, et_sta, et_ngb, et_nex, et_red, et_rec, et_suc, et_bad);
 //            printf("prepare:%.2f sec\tsearch:%.2f sec\tcheck:%.2f\n", lt_src, lt_pre, lt_pos);

}


//////////////////////////////////////////////////////////////////////
// Merge Log
//////////////////////////////////////////////////////////////////////

MergeLog::MergeLog()
{
    et_total = et_index = et_filter = et_search = et_path_index = et_simpaths = et_mtest = et_section = et_bases = et_align = et_update = 0.0;

    ct_allpaths = ct_simpaths = ct_bases = ct_align = ct_success = 0;
    ct_lgap = ct_egap = 0;
    ct_lgap_min = ct_egap_min = 1000000;
    ct_lgap_max = ct_lgap_sum = ct_egap_max = ct_egap_sum = 0;
	
    ct_bases_succ = ct_align_succ = 0;
    ct_align_succ_indel = 0;

    ct_align_fail_lgap = ct_align_fail_egap = ct_align_fail_short = ct_align_fail_score = 0;
    ct_align_succ_lfix = ct_align_succ_efix = ct_align_succ_lexp = ct_align_succ_eexp = 0;
}

void MergeLog::printSummary()
{
    printf("# total:%.4f\n", et_total);
    printf("# filtering:%.4f\n", et_filter);
    printf("# clustering:%.4f (path-search:%.4f, merge:%.4f)\n", et_search, et_simpaths, et_mtest);
    printf("- path-search:%.4f (path-index:%4f)\n", et_simpaths, et_path_index);
    printf("- merge:%.4f (merge-range:%.4f, simple-align:%.4f, alignments:%.4f)\n", et_mtest, et_section, et_bases, et_align);
    printf("- paths:%zu (all:%zu), simple-align:%zu (success:%zu), alignments:%zu (success:%zu [indel-align:%zu]), success:%zu\n", ct_simpaths, ct_allpaths, ct_bases, ct_bases_succ, ct_align, ct_align_succ, ct_align_succ_indel, ct_success );
    printf("- align failure with gaps (head-gap:%zu, tail-gap:%zu)\n", ct_align_fail_lgap, ct_align_fail_lgap);
    printf("- align failure with short length coverage:%zu\n", ct_align_fail_short);
    printf("- align failure with weak score:%zu\n", ct_align_fail_score);
    printf("- align success with fixed gaps (head:%zu, tail:%zu)\n", ct_align_succ_lfix, ct_align_succ_efix);
    printf("- align success with expansion (head:%zu, tail:%zu)\n", ct_align_succ_lexp, ct_align_succ_eexp);

    // printf("- head gaps (min:%zu, max:%zu, avg:%.2f)\n", ct_lgap_min, ct_lgap_max, (double)ct_lgap_sum/ct_lgap);
    // printf("- tail gaps (min:%zu, max:%zu, avg:%.2f)\n", ct_egap_min, ct_egap_max, (double)ct_egap_sum/ct_egap);
    printf("# update:%.4f\n", et_update);
    //printf("# index:%.4f\n", et_index);
}

//////////////////////////////////////////////////////////////////////
// Latch Log
//////////////////////////////////////////////////////////////////////

LatchLog::LatchLog()
{
    ljoin = 1;

    et_total = et_left = et_right = et_path = 0.0;
    et_path_filter = et_path_entry = et_path_sort = 0.0;
    et_reads = 0.0;
    et_latch = et_latch_check = 0.0;
    et_bases = 0.0;
    et_align = 0.0;
    et_connect = 0.0;
    et_section = 0.0;

    ct_path = ct_bases = ct_align = 0;
    ct_bases_succ = ct_align_succ = 0;
    ct_align_succ_indel = 0;
    ct_trial = ct_latch_succ = 0;
    
    ct_fail_pread = ct_fail_short = 0;
}

void LatchLog::setLJoin(bool l) 
{
    ljoin = l;
}

void LatchLog::setPESupport(bool p)
{
    pe_support = p;
}

void LatchLog::printSummary()
{
    printf("# total:%.4f, (left:%.4f, right:%.4f)\n", et_total, et_left, et_right );
    printf("# paths:%zu, time:%.4f\n", ct_path, et_path);
    if ( ljoin )
        printf("- kmer-filter:%.4f, entry-check:%.4f, sort:%.4f\n", et_path_filter, et_path_entry, et_path_sort );
    printf("# latch:%.4f (#trial:%zu, #success:%zu)\n", et_latch, ct_trial, ct_latch_succ);
    if ( pe_support ) 
        printf("- pair-reads:%.4f\n", et_reads );
    if ( ljoin ) 
        printf("- latchability:%.4f (section:%.4f), #base-comparisons:%zu (#success:%zu, time:%.4f) #aligns:%zu (#success:%zu [#indel-align:%zu], time:%.4f), connect:%.4f\n", et_section, et_latch_check, ct_bases, ct_bases_succ, et_bases, ct_align, ct_align_succ, ct_align_succ_indel, et_align, et_connect );
    else 
        printf("- connect:%.4f\n", et_connect );

    if ( pe_support ) 
        printf("extension fail (small pair reads:%zu, short paths:%zu)\n", ct_fail_pread, ct_fail_short);

}

BridgeLog::BridgeLog()
{
    pe_support = false;

    et_total = et_latchable = et_latch = 0.0;
    et_candidate = et_reads = et_drop = et_evidence = 0.0;
    et_update = 0.0;
    et_sort = 0.0;
    et_read_list = et_read_pivot = et_read_match = et_read_bad = et_read_trim = 0.0;
    ct_trial = ct_fail = ct_success = 0;
    ct_iter = 0;
    ct_candidates = ct_cand_max = 0;
    ct_cand_min = 1000000;
    ct_read_fail = ct_drop_fail = ct_evidence_fail = 0;
    ct_succ_first = ct_fail_first = 0;

	ct_read_all = ct_read_share = ct_read_align = ct_read_good = 0;
   
    ct_fail_pread = ct_fail_short = 0;
}

void BridgeLog::setPESupport( bool p )
{
    pe_support = p;
}

void BridgeLog::printSummary()
{
    printf("# total:%.4f, (find-paths:%.4f [sort:%.4f], latchablity:%.4f, latch:%.4f, update:%.4f)\n", et_total, et_candidate, et_sort, et_latchable, et_latch, et_update );
    printf("- latchability (find-reads:%.4f, trim-reads:%.4f, evidence:%.4f)\n", et_reads, et_drop, et_evidence);
    printf("- find reads (pivot-map:%.4f, match-map:%.4f, bad-reads:%.4f, trim-reads:%.4f)\n", et_read_pivot, et_read_match, et_read_bad, et_read_trim );
    if ( pe_support ) 
        printf("- pair reads:%.4f\n", et_read_pair);
    printf("# bridigng reads (all:%zu, shared:%zu, aligned:%zu, good:%zu)\n", ct_read_all, ct_read_share, ct_read_align, ct_read_good);
    printf("# candidate paths:%zu (max:%zu, min:%zu, avg:%.2f)\n", ct_candidates, ct_cand_max, ct_cand_min, (ct_candidates/(double)ct_iter));
    //printf("- failure (small-reads:%zu, trimmed-reads:%zu, poor-evidence:%zu)\n", ct_read_fail, ct_drop_fail, ct_evidence_fail );
    printf("# iteration:%zu, (success:%zu, fail:%zu)\n", ct_trial, ct_success, ct_fail );
    printf("# top candidate only (success:%zu, fail:%zu)\n", ct_succ_first, ct_fail_first );
    printf("- failure (small-reads:%zu, poor-evidence:%zu)\n", ct_read_fail, ct_evidence_fail );
    printf("- failure (pair-reads:%zu, short-paths:%zu)\n", ct_fail_pread, ct_fail_short );
}

//////////////////////////////////////////////////////////////////////
// Recruit Log
//////////////////////////////////////////////////////////////////////
RecruitLog::RecruitLog()
{
    by_anchor = true;

    ct_sub = ct_aln = 0;
    ct_cmp = ct_good = 0;
    ct_can = 0;
    et_sub = et_aln = et_pos = 0.0;
    et_anc = 0.0;
    et_can = 0.0;
    sub_good = aln_good = 0;

    et_can_all = et_can_good = 0.0;
    et_sbjct_find = et_sbjct_get = 0.0;
}

void RecruitLog::print()
{
    printf("# Candidates: %zu\n", ct_can);
    if ( by_anchor ) {
        printf("# Comparisons: %zu\n", ct_cmp );
        printf("# Recruited: %zu\n", ct_good);
        printf("# Search:%.2f (paths:%.2f, trim:%.2f), anchoring:%.2f, sbjct-find:%.2f sbjct-get:%.2f\n", et_can, et_can_all, et_can_good, et_anc, et_sbjct_find, et_sbjct_get);
    } else {
        printf("# Comparisons: %zu (substring:%zu, alignment:%zu)\n", 
               ct_sub+ct_aln, ct_sub, ct_aln );
        printf("# Recruited: %zu (substring:%zu, alignment:%zu)\n",
               sub_good+aln_good, sub_good, aln_good);
        
        printf("# Search:%.2f, substring:%.2f, range-search:%.2f, alignment:%.2f\n",
               et_can, et_sub, et_pos, et_aln);
    }

}

//////////////////////////////////////////////////////////////////////
// Anchor Log
//////////////////////////////////////////////////////////////////////

AnchorLog::AnchorLog()
{
    et_init = et_find = et_range = et_map = et_substr = et_score = 0.0;
    et_init_read = et_init_path = et_init_insert = 0.0;
    ct_all = ct_valid = ct_exist = ct_pass = 0;
}

void AnchorLog::add( AnchorLog &other )
{
    et_init += other.et_init;
    et_init_insert += other.et_init_insert;
    et_init_read += other.et_init_read;
    et_init_path += other.et_init_path;
    et_find += other.et_find;
    et_range += other.et_range;
    et_map += other.et_map;
    et_substr += other.et_substr;
    et_score += other.et_score;

    ct_all += other.ct_all;
    ct_valid += other.ct_valid;
    ct_exist += other.ct_exist;
    ct_pass += other.ct_pass;
}

void AnchorLog::print()
{
    printf("# Anchoring summary:\n");
    printf("# all:%zu, valid-range:%zu, dup-range:%zu, good-filter:%zu\n", ct_all, ct_valid, ct_exist, ct_pass);
    //printf("# init:%.4f (read:%.4f, path:%.4f, insert:%.4f)\n", et_init, et_init_read, et_init_path, et_init_insert);
    printf("# init:%.4f\n", et_init);
    printf("# find:%.4f\n", et_find);
    printf("- range:%.4f, map:%.4f, substr:%.4f, score:%.4f\n", et_range, et_map, et_substr, et_score);
}
