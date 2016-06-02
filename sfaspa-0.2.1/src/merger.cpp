#include "merger.h"

Merger::Merger()
{
    kmer_size = Param::merge_filter_kmer;

    status = false;
}

Merger::~Merger()
{
    clear();
}

void Merger::init( PathEntryMap* paths )
{
    double t0 = mytime();

    path_entries = paths;

    PathId pivot = getMaxPathId();
    if ( Param::verbose ) std::cout << "Max pid:" << pivot << "\n";
    if ( Param::verbose ) std::cout << "Size: " << path_entries->size() << "\n";

    PathEntryMap::const_iterator it;
    for ( it = path_entries->begin(); it != path_entries->end(); ++it ) {
        //-----------------------
        // Check numeric ID limit
        //-----------------------
        assert( pivot < std::numeric_limits<PathId>::max() );
        ++pivot;
        
        //-------------------------------------------------
        // Create an connect path object from a path entry
        //-------------------------------------------------
        PathId pid = it->first;          // existing Path ID
        std::string str = it->second.seq; // path sequence
        AlignSummary summary;
        summary.self(str.size());

        Cluster c;
        c.sequence = str;
        c.add(pivot, summary);
        c.lstop = it->second.lstop;
        c.rstop = it->second.rstop;
        c.ltrim = it->second.ltrim;
        c.rtrim = it->second.rtrim;
        c.lbase = 0;
        c.rbase = 0;

        //-----------------------------
        // Save the just created object
        //-----------------------------
        clusters.insert( std::pair<PathId, Cluster>(pivot, c) );

        //--------------------------------------------
        // Path ID mapping between of old/new path ID
        //--------------------------------------------
        id_map.insert( std::pair<PathId, PathId>( pivot, pid ) );
    }
    
    if ( Param::verbose >= 1 ) 
        std::cout << "Clustering entries initialized:" << mytime()-t0 << " sec\n";
}


void Merger::clear()
{
    path_entries = NULL;

    kmer_paths.clear();
    path_kmers.clear();
    merged_paths.clear();
    clusters.clear();

    status = false;

    if ( Param::verbose ) std::cout << "Merger cleaned\n";

}

void Merger::release()
{
    kmer_paths.clear();
    path_kmers.clear();
}

void Merger::run()
{
    build();
    cluster();
    update();
    dropMergedPaths();
    writeClusters();
    release();

    status = true;
}

void Merger::build()
{
    PathEntryMap::const_iterator it;
    //for ( it = path_entries.begin(); it != path_entries.end(); ++it ) {
    for ( auto it = clusters.begin(); it != clusters.end(); ++it ) {
        PathId pid = it->first;
        //if ( merged_paths.find(pid) != merged_paths.end() ) continue;
        std::string str = it->second.sequence;
        size_t len = str.size();
        assert(len>0);

        KmerArray kids = biostr::getKmers(str, kmer_size);
        path_kmers.add(pid, kids);

        KmerPosMap poss;
        for ( size_t j = 0; j < kids.size(); j++ )
            poss[kids[j]].push_back(j);
        path_kposs.add(pid,poss);
    }
}

void Merger::cluster()
{
    double t0    = mytime();
    
    PathLengthMap plen_map = getPivotLengths();
    size_t npath = plen_map.size();

    progress = Progress( 1.0, 0, npath, mytime() );

    PathLengthMap::reverse_iterator it;
    for ( it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId query_pid = it->second;

        assert( merged_paths.find(query_pid) == merged_paths.end() );

        double prev = progress.ratio;
        log.et_total = (mytime()-t0);

        /* Print progress status at every percent of path proceeded */
        progress.count++;        
        progress.showProgress();
        if ( progress.ratio > prev ) {
            if ( Param::verbose >= 1 ) 
                log.printSummary();
        }
        if ( Param::verbose ) {
            //std::string query_str = path_entries[query_pid].seq;
            std::string query_str = getPathSequence(query_pid);
            std::cout << "\nQuery:" << query_pid 
                      << "\tLength:" << query_str.size() << "\n"
                      << query_str << "\n";
        }
        
        PathId sbjct_pid;
        AlignSummary aln_summary;
        
        //std::string query = path_entries[query_pid].seq;
        std::string query = getPathSequence(query_pid);

        if ( ! Param::merge_fixed_kmer ) {
            int ksize = getNewFilterKmer( query.size() );

            if ( ksize < kmer_size ) {
                if ( ksize < 2 ) {
                    if ( Param::verbose >= 1 ) 
                        std::cout << "Too small k-mer:" << ksize << "\tsequence length:" << query.size() << "\n";
                    break;
                }
                if ( Param::verbose >= 1 ) {
                    std::cout << "New filter k-mer:" << ksize << "\tsequence length:" << query.size() << "\n";
                }
                kmer_size = ksize;
                
                rebuildIndex();
            }
        }

        //double tic = mytime();
        //KmerArray  query_kmers = biostr::getKmers(query, Param::merge_filter_kmer);
        //KmerArray  query_kmers = biostr::getKmers(query, kmer_size);
        const KmerArray*  query_kmers = path_kmers.getPointer(query_pid);
        //log.et_kmers += (mytime()-tic );
        assert( query_kmers->size() > 0 );
        
        double tic = mytime();
        bool found = findCluster( query_kmers, query_pid, sbjct_pid, aln_summary );
        log.et_search += (mytime()-tic);

        tic = mytime();
        if ( found ) {
            merged_paths.insert( query_pid );
            //clusters[sbjct_pid].push_back( PathAlignPair(query_pid, aln_summary) );
            clusters[sbjct_pid].add(query_pid, aln_summary);
            if ( Param::verbose ) std::cout << "success:"  << mytime()-t0 << " sec\tsbjct:" << sbjct_pid << "\tsize:" << clusters[sbjct_pid].members.size() << "\n";
        }
        else {
            // aln_summary.init(); // reset alignment 
            // aln_summary.self(query.size()); // self alignment

            // //clusters[query_pid].push_back( PathAlignPair(query_pid, aln_summary) );
            // clusters[query_pid].add(query_pid, aln_summary);

            // Cluster has been already created during initialization
            kmer_paths.add( *query_kmers, query_pid );
        }
        log.et_update += (mytime()-tic);
    }

    if ( Param::verbose ) std::cout << "# Merged paths:" << merged_paths.size() << "\n";

    if ( Param::verbose >= 1 ) {
        std::cout << "# Candidate:" << log.ct_simpaths << "\n";
        std::cout << "# Base comparison:" << log.ct_bases << "\n";
        std::cout << "# Pairwise Alignment:" << log.ct_align << "\n";
        std::cout << "# Success:" << log.ct_success << "\n";
        // std::cout << "# Head gap fail:" << log.ct_align_fail_lgap << "\n";
        // std::cout << "# Tail gap fail:" << log.ct_align_fail_egap << "\n";
        std::cout << "# Sum align time:" << log.et_align << "\n";
        std::cout << "# Avg align time:" << log.et_align / log.ct_align << "\n";
    }
    std::cout << "# Merged:" << log.ct_success << "\n";
    //std::cout << "# Clusters:" << clusters.size() << "\n";
    std::cout << "# Clusters:" << npath - merged_paths.size() << "\n";
}

int Merger::getNewFilterKmer( size_t str_len )
{
    double t0 = mytime();
    int ksize = filter::getFilterKmer( str_len, Param::merge_shared_nkmer, 1.0-Param::merge_filter_score );
    
    // if ( Param::use_preset_filter ) {
    //     ksize = kmer_size;
    //     int nkmer;
    //     double score;
    //     while (ksize>=2) {
    //         score = prefilter.getFilter( ksize );
    //         nkmer = filter::minSameKmerCount( str_len, ksize, 1-score);
    //         if ( nkmer < Param::merge_shared_nkmer ) ksize--;
    //         else break;
    //     }
    // }
    log.et_filter += (mytime()-t0);
    return ksize;
}

void Merger::findSimilarPaths( //PathFreqMap &path_freq,
                              int min_filter,
                              std::multimap<size_t, PathId> &freq_map,
                              PathId &query_pid,
                              std::string &query,
                              const KmerArray* query_kmers )
{
    double t0 = mytime();        
    

    PathFreqMap path_freq;
    kmer_paths.findPaths( path_freq, query_kmers );    
    log.et_path_index += (mytime()-t0);
    log.ct_allpaths += path_freq.size();

    //std::multimap<size_t, PathId> freq_map;
    //std::multimap<size_t, PathId>::reverse_iterator it;
    for ( PathFreqMap::iterator pt = path_freq.begin(); pt != path_freq.end(); ++pt ) {
        if ( pt->first == query_pid ) continue;
        if ( (int)pt->second < min_filter ) continue;
        freq_map.insert( std::pair<size_t, PathId>(pt->second, pt->first) );
    }


    if ( Param::verbose ) 
        std::cout << "#Paths:" << freq_map.size() 
                  << "\tScan paths time:" << mytime() -t0 << "\n";
}

bool Merger::mergible( int min_filter,
                       std::string &sbjct,
                       std::string &query, 
                       const KmerArray* skmers,
                       const KmerArray* qkmers,
                       const KmerPosMap *sposs,
                       const KmerPosMap *qposs,
                       AlignSummary &summary )
{
    //-----------------------------
    // Initialize alignment summary
    //-----------------------------
    summary.init();

    //---------------------------------------------
    // Find mergible region between query and sbjct
    //---------------------------------------------
    double tic = mytime();
    MergeSectioner sect(&query, &sbjct, qkmers, skmers, qposs, sposs);
    sect.find();
    Section s = sect.getSection();
    assert( s.qbeg == 0 && s.qend == (int)query.size()-1 ); // full query region
    log.et_section += (mytime()-tic);

    //---------------------------
    // Do base by base comparison
    //---------------------------
    if ( Param::align_base_first ) 
        if ( doBaseAlign( summary, query, sbjct, s ) ) return true;
    
    //--------------------------------------------------------
    // Make sbjct region a bit wider before pairwise alignment
    //--------------------------------------------------------
    sect.pad();
    s = sect.getSection();
    
    //---------------------------
    // Perform pairwise alignment
    //---------------------------
    return doPairAlign( summary, query, sbjct, s );
}

bool Merger::doBaseAlign( AlignSummary &summary,
                          std::string &query,
                          std::string &sbjct,
                          Section &s )
{
    //---------------------------------------------------
    // Must be same length between query and sbjct region
    //---------------------------------------------------
    if ( s.send-s.sbeg != s.qend-s.qbeg ) return false;
    
    int qs = s.qbeg, ss = s.sbeg;
    int l = s.qend-s.qbeg+1;

    assert(qs==0);
    assert(ss>=0);
    assert(l>0);

    double tic = mytime();
    double align_score = scoring::positiveRate( &query, &sbjct, qs, ss, l, BLOSUM62 );
    log.et_bases += (mytime()-tic);
    log.ct_bases++;
    
    if ( align_score >= Param::merge_score ) {
        if ( Param::verbose ) {
            printf("score:%.2f\n", align_score);
            printf("expand:\tqbeg:%d sbeg:%d qend:%d send:%d\n", s.qbeg,s.sbeg, s.qend,s.send);
        }
        
        log.ct_success++;
        log.ct_bases_succ++;
        
        summary.score = align_score;
        if ( s.sbeg > 0 ) 
            summary.lgap.second = s.sbeg;
        
        if ( Param::verbose ) {
            std::cout << "Simple align success\n";
            summary.print(std::cout);
        }
        return true;
    }
    return false;
}

bool Merger::doPairAlign( AlignSummary &summary,
                          std::string &query,
                          std::string &sbjct,
                          Section &s )
{
    //----------------------------------------------------------------
    // Clean up alignment summary here, again, because it was modified
    // during base/base comparison
    //----------------------------------------------------------------
    summary.init();
    
    if ( Param::verbose ) {
        std::cout << "sbeg:" << s.sbeg << "\tsend:" << s.send << "\n";
        std::cout << "qbeg:" << s.qbeg << "\tqend:" << s.qend << "\n";
    }

    std::string sbjct_sub = sbjct.substr(s.sbeg, s.send-s.sbeg+1);
    //std::string query_sub = query.substr(s.qbeg, s.qend-s.qbeg+1);
    
    
    if ( Param::verbose ) {
        std::cout << "Sbjct sub:" << sbjct_sub << "\n";
        //std::cout << "Query sub:" << query_sub << "\n";
    }

    double tic = mytime();

    int band_size = int ( Param::band_ratio * sbjct_sub.size() );
    SemiGlobalAlign aln = ! Param::banded_align ? 
        SemiGlobalAlign(sbjct_sub, query, ANCHOR_CENTER, Param::gap_ext, Param::gap_open) :
        SemiGlobalAlign(sbjct_sub, query, ANCHOR_CENTER, Param::gap_ext, Param::gap_open, -1*band_size, band_size) ;

    log.et_align += (mytime()-tic);
    log.ct_align ++;

    summary = aln.getSummary();
    
    
    if ( Param::verbose ) {
        aln.printAlignment(std::cout);
        summary.print(std::cout);
        summary.printGap(std::cout);
    }

    //-------------------------------------------------------------------
    // Invalid alignmenht (both sbjct and query has non-zero leading gaps
    //-------------------------------------------------------------------
    if ( summary.lgap.first > 0 && summary.lgap.second > 0 ) {
        if ( Param::verbose ) std::cout << "Leading gaps in both paths\n";
        return false;
    }


    //--------------------------------------------------------------------
    // Invalid alignmenht (both sbjct and query has non-zero trailing gaps
    //--------------------------------------------------------------------
    if ( summary.egap.first > 0 && summary.egap.second > 0 ) {
        if ( Param::verbose ) std::cout << "Trailing gaps in both paths\n";
        return false;
    }

    //--------------------
    // Get alignment score
    //--------------------
    double score = !Param::identity_flag ? 
        summary.posrate :
        (double)summary.match/summary.length;

    //-----------
    // Weak score
    //-----------
    if ( score < Param::merge_score ) {
        if ( Param::verbose ) 
            std::cout << "Weak score:" << score << "\n";
        log.ct_align_fail_score++;
        return false;
    }
    
    //-----------------------
    // Short alignment region
    //-----------------------
    double ratio = (double)summary.length / query.size();
    if ( ratio < Param::merge_short_ratio ) {
        if ( Param::verbose ) 
            std::cout << "Short aligned ratio:" << ratio << "\n";
        log.ct_align_fail_short++;
        return false;
    }

    //--------------------------------------------------------------
    // sbjct sequence head offset: sbjct region start - leading gaps
    // xxxxxoooooooooooooooxxxxx sbjct
    //    --ooooooooooooooo      sbjct sub
    //    xxooooooooooooooo      query
    // sbeg = 5, lgap = 2
    // hoff = 3
    //--------------------------------------------------------------
    int hoff = s.sbeg - summary.lgap.first;

    //--------------------------------------------------------------
    // sbjct sequence tail offset
    // xxxxxoooooooooooooooxxxxx sbjct
    //      ooooooooooooooo--    sbjct sub
    //      oooooooooooooooxx    query
    // size:25, s.send:19, egap = 2
    // toff =  (25-20)-2 = 3
    //--------------------------------------------------------------
    int toff = ( (int)sbjct.size() - (s.send+1) ) - summary.egap.first;

    //================================================================
    // Case: leading gaps in sbjct sequence
    //================================================================
    if ( summary.lgap.first > 0 ) { 
        assert( summary.lgap.second == 0 );

        //---------------------------------
        //  xoooooooooooooooxxxxx sbjct
        //---ooooooooooooooo      sbjct sub
        //xxxooooooooooooooo      query
        // hoff = -2
        //---------------------------------
        if ( hoff < 0 ) {
            int diff = summary.lgap.first-s.sbeg;
            if ( diff < (int)log.ct_lgap_min ) log.ct_lgap_min = diff;
            if ( diff > (int)log.ct_lgap_max ) log.ct_lgap_max = diff;                
            log.ct_lgap_sum += diff;
            log.ct_lgap++;
            
            //-------------------------------------------
            // Do not allow sbjct sequence to expand left
            //-------------------------------------------
            if ( ! Param::extend_cluster ) {
                if ( Param::verbose ) std::cout << "Leading gap -> skip\n";
                log.ct_align_fail_lgap++;
                return false; 
            } 
            //--------------------------------------
            // Allow sbjct sequence expansion (left)
            //--------------------------------------
            else {
                if ( Param::verbose ) std::cout << "soff:" << hoff << "\n";
                //-------------------------------------
                // Adjust gaps & positions
                //------------------------
                //    xxxooooooooooxxxxxxxxx
                // -----xooooooooooxxxxxxxxx sbjct sub
                // xxxxxxoooo-ooooo          query
                // hoff = -3
                // s.sbeg = 2
                // sbjct lgap = 5 -> 3
                // sbjct indel : 5 -> 7
                //------------------------------------
                summary.lgap.first = abs(hoff); // sbjct leading gap
                //------------------------------------------------------------
                // Since query region is the entire sequence range, 
                // it does not need to update INDES postions in query.
                // But, it does need to update INDELs in sbjct by moving right.
                //------------------------------------------------------------
                if ( s.sbeg > 0 )
                    summary.shift(s.sbeg, true); // adjust sbjct postions of INDELs
                log.ct_align_succ_lexp++;
                if ( Param::verbose ) 
                    std::cout << "Good aligned ratio:" << ratio << "\n";
            }
        } 
        //------------------------------------
        // Case: Leading gap correction
        // head off >= 0
        //------------------------------------
        else {
            if ( Param::verbose ) std::cout << "Leading gap -> fixable\n";

            //------------------------------------
            // xxxxxxooooooooooxxxxxxxxx
            //  ----xooooooooooxxxxxxxxx sbjct sub
            //  xxxxxoooo-ooooo          query
            // hoff = 1
            // lgap = 4
            // sbeg = 5
            // sbjct lgap => 0
            // query lgap => 1 (hoff)
            //------------------------------------
            summary.lgap.first = 0;
            summary.lgap.second = hoff;
            
            //--------------
            // Adjsut INDELs
            //-------------- 
            if ( s.sbeg > 0 ) summary.shift(s.sbeg, true); 

            log.ct_align_succ_lfix++;
            if ( Param::verbose ) std::cout << "New sbeg:" << s.sbeg << "\n";
        }
    }
    
    //================================================================
    // Case: NO leading gaps in sbjct sequence
    //================================================================
    else {
        // xxxxxooooooooooxxxxx sbjct
        //    xxoooooooooo      sbjct sub
        //    --oooooooooo      query
        // s.sbeg = 3
        //----------------------------------------
        // adjust leading gaps and INDEL positions
        //----------------------------------------
        if ( s.sbeg > 0 ) {
            summary.lgap.second += s.sbeg;
            summary.shift(s.sbeg, true);
        }
    }

    //================================================================
    // Case: trailing gaps in sbjct sequence
    //================================================================
    if ( summary.egap.first > 0 ) { 
        assert( summary.egap.second == 0 );

        //--------------------------------------------------------------
        // sbjct sequence tail offset
        // xxxxxooooooooooooooox   sbjct
        //      ooooooooooooooo-   sbjct sub
        //      oooooooooooooooxxx query
        // toff = -2
        //--------------------------------------------------------------
        if ( toff < 0 ) {
            int diff = summary.egap.first - ((int)sbjct.size() - (s.send+1)) ;
            if ( diff < (int)log.ct_egap_min ) log.ct_egap_min = diff;
            if ( diff > (int)log.ct_egap_max ) log.ct_egap_max = diff;                
            log.ct_egap_sum += diff;
            log.ct_egap++;
            
            //--------------------------------------------
            // Do not allow sbjct sequence to expand right
            //--------------------------------------------
            if ( ! Param::extend_cluster ) {
                if ( Param::verbose ) std::cout << "Trailing gap -> skip\n";
                log.ct_align_fail_egap++;
                return false; 
            }
            
            //---------------------------------------
            // Allow sbjct sequence expansion (right)
            //---------------------------------------
            else {
                //--------------------------
                // Adjust gaps.
                // No need to update INDELs
                //-------------------------
                summary.egap.first = abs(toff); // sbjct trailing gap
                log.ct_align_succ_eexp++;
                if ( Param::verbose ) 
                    std::cout << "Good aligned ratio:" << ratio << "\n";
            }
        } 
        //------------------------------------
        // Case: Trailing gap correction 
        // tail toff >= 0
        //------------------------------------
        else {
            //-------------------------------
            // xxxxxoooooooooooxxxx sbjct
            //      ooooooooooo---  sbjct sub
            //      oooooooooooxxx  query
            // toff = 1
            // egap = 3
            // sbjct egap => 0
            // query lgap => 1 (hoff)
            //-------------------------------

            if ( Param::verbose ) std::cout << "Trailing gap -> fixable\n";
            summary.egap.first  = 0; // sbjct trailing gap
            summary.egap.second = toff;
            log.ct_align_succ_efix++;
        }
    }
    //================================================================
    // Case: NO trailing gaps in sbjct sequence
    //================================================================
    else {
        int toff = ( (int)sbjct.size() - (s.send+1) );
        summary.egap.second += toff;
    }

    log.ct_success++;
    log.ct_align_succ++;
    

    if ( Param::verbose ) {
        std::cout << "Alignment summary updated\n";
        summary.print(std::cout);
        summary.printGap(std::cout);
    }
    
    if ( summary.ilist.size() || summary.dlist.size() ) log.ct_align_succ_indel++;
    
    return true;
}

bool Merger::findCluster( const KmerArray* query_kmers,
                          PathId &query_pid, 
                          PathId &sbjct_pid, 
                          AlignSummary &summary )
{
    //std::string query = path_entries[query_pid].seq;
    std::string query = getPathSequence(query_pid);
    
    //int min_filter = filter::minSameKmerCount( query.size(), Param::merge_filter_kmer, 1-Param::merge_filter_score );
    int min_filter = filter::minSameKmerCount( query.size(), kmer_size, 1-Param::merge_filter_score );
    //if ( min_filter < 1 ) min_filter = 1;
    if (Param::verbose) std::cout << "Min k-mer:" << min_filter << "\n";
    if ( min_filter < Param::merge_shared_nkmer ) {
        if ( Param::verbose ) std::cout << "Small k-mer match\n";
        return false;
    }
    
    double tic = mytime();
    std::multimap<size_t, PathId> freq_map;
    findSimilarPaths( min_filter, freq_map, query_pid, query, query_kmers );
    log.et_simpaths += (mytime()-tic);
    log.ct_simpaths += (freq_map.size());
    if ( freq_map.size() ==  0 ) return false;

    //ount_pat += freq_map.size();


    const KmerPosMap* qposs = path_kposs.getPointer(query_pid);

    // //KmerArray qkmers = biostr::getKmers(query, kmer_size);
    // KmerPossMap qposs;
    // for ( size_t i = 0; i < query_kmers.size(); i++ ) 
    //     qposs[query_kmers[i]].push_back(i);

    size_t niter = 0;
    std::multimap<size_t, PathId>::reverse_iterator it;
    for ( it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
        sbjct_pid = it->second;
        //std::string sbjct = (*path_entris)[sbjct_pid].seq;
        std::string sbjct = getPathSequence(sbjct_pid);
        assert(sbjct.size()>0);
        
        if ( Param::verbose ) 
            std::cout << "Iteration:" << ++niter 
                      << "\tSbjct:" << sbjct_pid 
                      << "\tLength:" << sbjct.size() 
                      << "\t# Shared kmers:" << it->first << "\n"
                      << sbjct << "\n";

        //KmerArray skmers = biostr::getKmers(sbjct, kmer_size);        
        const KmerArray* sbjct_kmers = path_kmers.getPointer(sbjct_pid);
        
        assert(sbjct_kmers->size() > 0 );
        const KmerPosMap *sposs = path_kposs.getPointer(sbjct_pid);
        // KmerPossMap sposs;
        // for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) 
        //     sposs[sbjct_kmers[i]].push_back(i);
        
        double tic = mytime();
        bool good = mergible( min_filter, sbjct, query, sbjct_kmers, query_kmers, sposs, qposs, summary );
        log.et_mtest += (mytime()-tic);
        if ( good ) return true;
    }
    return false;
}

void Merger::update()
{
    for ( auto it = clusters.begin(); it != clusters.end(); ++it ) {
        PathId pivot = it->first;
        if ( merged_paths.find(pivot) != merged_paths.end() ) continue;
        
        size_t max_lgap = 0, max_egap = 0;
        PathId max_lpid, max_epid;
        
        it->second.findMaxGap( max_lgap, max_lpid, 1 );
        it->second.findMaxGap( max_egap, max_epid, 0 );

        if ( max_lgap > 0 ) {
            // if ( clusters.find(max_lpid) == clusters.end() ) {
            //     std::cout << "Somthing wrong:max-lpid:" << max_lpid << "\tmax_lgap:" << max_lgap << "\n";
            // }
            assert( clusters.find(max_lpid) != clusters.end() );
            std::string str = clusters[max_lpid].sequence;
            assert( str.size() > max_lgap );
            it->second.sequence = str.substr(0, max_lgap) + it->second.sequence;
            it->second.lbase = max_lgap;
        }

        if ( max_egap > 0 ) {
            // if ( clusters.find(max_epid) == clusters.end() ) {
            //     std::cout << "Somthing wrong:max-epid:" << max_epid << "\tmax_egap:" << max_egap << "\n";
            // }

            assert( clusters.find(max_epid) != clusters.end() );
            std::string str = clusters[max_epid].sequence;
            assert( str.size() > max_egap );
            it->second.sequence += str.substr(str.size()-max_egap);
            it->second.rbase = max_egap;
        }

        if ( max_lgap ) 
            it->second.updatePositions( max_lgap );
    }
}

void Merger::writeClusters()
{
    if ( ! Param::output_all ) return;

    //if ( Param::extend_first ) assert( Param::recluster_flag == false );

    std::string file;
    char htype, mtype;
    if ( !recall ) {
        file = Param::out_dir + "/cluster.fasta";
        htype = 'c';
        if ( Param::extend_first ) 
            mtype = ( Param::read_bridge_extend ) ? 'b' : 's';        
        else 
            mtype = 'p';
    } 
    // else {
    //     file = Param::out_dir + "/recluster.fasta" ;
    //     htype = 'C';
    //     mtype = ( Param::read_bridge_extend ) ? 'b' : 's';        
    // }

    std::fstream out;
    fio::openFile( out, file.c_str(), std::ios::out );
    

    PathLengthMap plen_map = getPivotLengths();
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;

        auto jt = clusters.find(pid);
        assert( jt != clusters.end() );

        std::string seq = jt->second.sequence;
        out << ">" << htype << pid
            << " len:" << seq.size() 
            << " stop:" << jt->second.lstop << "/" << jt->second.rstop
            << " trim:" << jt->second.ltrim << "/" << jt->second.rtrim
            << " expand:" << jt->second.lbase << "/" << jt->second.rbase
            << " mem:" << clusters[pid].members.size() << " ";
        for ( auto mem : clusters[pid].members ) {
            out << mtype;
            PathId oid = id_map[mem];
            out << oid << ";";
        }
        
        out << "\n";
        out << seq << "\n";
    }
}

// PathLengthMap Merger::getPathLengths()
// {
//     PathLengthMap lmap;

//     PathEntryMap::const_iterator it;
//     for ( it = path_entries.begin(); it != path_entries.end(); ++it ) {
//         PathId pid = it->first;
//         if ( merged_paths.find(pid) != merged_paths.end() ) continue;
//         std::string str = it->second.seq;
//         size_t len = str.size();
//         assert(len>0);
//         lmap.insert( std::pair<size_t, PathId>( len, pid ) );
//     }

//     return lmap;
// }

PathLengthMap Merger::getPivotLengths()
{
    PathLengthMap lmap;

    for ( auto it = clusters.begin(); it != clusters.end(); ++it ) {
        PathId pid = it->first;
        if ( merged_paths.find(pid) != merged_paths.end() ) continue;
        assert( clusters.find(pid) != clusters.end() );
        std::string str = it->second.sequence;
        size_t len = str.size();
        assert(len>0);
        lmap.insert( std::pair<size_t, PathId>( len, pid ) );
    }

    return lmap;
}

void Merger::print()
{
    int count = 0;
    PathLengthMap plen_map = getPivotLengths();
    for ( PathLengthMap::reverse_iterator it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;
        assert(pid>=0);
        //int pos = (int)pid;
        assert( merged_paths.find(pid) == merged_paths.end() );
        assert( clusters.find(pid) != clusters.end() );
        
        std::cout << ++count << "\tRepresentative:" << pid << "\n";
        std::string ref = clusters[pid].sequence;
        std::cout << ref << "\n";
        std::cout << "# Members:" << clusters[pid].members.size() << "\n";
        
        PathIdList members = clusters[pid].members;
        AlignSummaryList summarys = clusters[pid].summarys;
        
        PathIdList::iterator ct;
        AlignSummaryList::iterator st;
        for ( ct = members.begin(), st = summarys.begin(); 
              ct != members.end() && st != summarys.end();
              ++ct, ++st ) {

            //for ( auto jt = clusters[pid].begin(); jt != clusters[pid].end(); ++jt ) {
            //PathId mem = jt->first;
            PathId mem = *ct;
            int pos = st->lgap.second;
            //int pos = jt->second.lgap.second;
            std::string mstr = getPathSequence(mem);
            std::cout << mem << "\t" << pos << "\t" << mstr << "\n";
        }
        std::cout << "----------------------------------------------------------------------\n";
    }
}

void Merger::dump( std::string filename )
{
    std::fstream out;
    fio::openFile( out, filename.c_str(), std::ios::out | std::ios::binary );

    int count = 0;
    PathLengthMap plen_map = getPivotLengths();
    size_t size = plen_map.size();
    out.write((char*)&size, sizeof(size_t));
    PathLengthMap::reverse_iterator it;
    for ( it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        PathId pid = it->second;
        assert(pid>=0);

        assert( merged_paths.find(pid) == merged_paths.end() );
        assert( clusters.find(pid) != clusters.end() );

        out.write((char*)&pid, sizeof(PathId));

        // size_t len = sequence.size();
        // const char *cstr = sequence.c_str();
        // out.write((char*)&len, sizeof(size_t));
        // out.write(cstr, len);

        clusters[pid].dump(out);
        // size_t ctmem = clusters[pid].members.size();
        // out.write((char*)&ctmem, sizeof(size_t));
        // for ( auto item : clusters ) 
        //     item.dump(out);

        // PathAlignPairList::iterator jt;
        // for ( jt = clusters[pid].begin(); jt != clusters[pid].end(); ++jt ) {
        //     PathId mem = jt->first;
        //     AlignSummary aln = jt->second;
        //     out.write((char*)&mem, sizeof(PathId));
        //     aln.dump(out);
        // }
        
        //path_entries[pid].dump(out);

        ++count;
    }
    
    size = merged_paths.size();
    out.write((char*)&size, sizeof(size_t));
    for ( PathIdSet::iterator it = merged_paths.begin(); it != merged_paths.end(); ++it )
        out.write((char*)&(*it), sizeof(PathId));
    

    size = id_map.size();
    out.write((char*)&size, sizeof(size_t));
    for ( PathIdMap::iterator it = id_map.begin(); it != id_map.end(); ++it ) {
        out.write((char*)&(it->first),  sizeof(PathId));
        out.write((char*)&(it->second), sizeof(PathId));
    }

    if ( Param::verbose ) std::cout << "#clusters dumped:" << count << "\n";
}

void Merger::load( std::string filename )
{
    std::fstream in;
    fio::openFile( in, filename.c_str(), std::ios::in | std::ios::binary );

    int count = 0;
    size_t size;
    in.read((char*)&size, sizeof(size_t));
    for ( size_t i = 0; i < size; i++ ) {
        PathId pid;
        in.read((char*)&pid, sizeof(PathId));

        // size_t len;
        // in.read((char*)&len, sizeof(size_t));
        // assert(len>0);
        // char *seq = new char[len+1];
        // in.read(seq, len);
        // seq[len]='\0';
        // sequence = std::string(seq);
        // delete[] seq;
        
       //  size_t ctmem;
       //  in.read((char*)&ctmem, sizeof(size_t));
       //  clusters.insert(std::pair<PathId, PathAlignPairList>(pid, PathAlignPairList()));
       //  for ( size_t j = 0; j < ctmem; j++ ) {
       //      Cluster c;
       //      c.load(
       //      PathId mem;
       //      AlignSummary aln;
       //      in.read((char*)&mem, sizeof(PathId));
       //      aln.load(in);
       //      clusters[pid].push_back( PathAlignPair( mem, aln ) );
       // }

        Cluster c;
        c.load(in);
        clusters.insert( std::pair<PathId, Cluster>( pid, c) );

        // PathEntry entry;
        // entry.load(in);
        // path_entries.insert( std::pair<PathId, PathEntry>( pid, entry ) );
        
        ++count;
    }

    PathId pid;
    in.read((char*)&size, sizeof(size_t));
    for ( size_t i = 0; i < size; i++ ) {
        in.read((char*)&pid, sizeof(PathId));
        merged_paths.insert(pid);
    }
        
    PathId nid, oid;
    in.read((char*)&size, sizeof(size_t));
    for ( size_t i = 0; i < size; i++ ) {
        in.read((char*)&nid, sizeof(PathId));
        in.read((char*)&oid, sizeof(PathId));
        id_map.insert( std::pair<PathId, PathId>(nid, oid) );
    }

    printElapsed( INIT_TIME, mytime(), "Clusters loaded" );

    status = true;
}

void Merger::dropMergedPaths()
{
    for ( Clusters::iterator it = clusters.begin(); it != clusters.end(); ) 
        merged_paths.find(it->first) == merged_paths.end() ?
            ++it : clusters.erase(it++);
}

// //void Merger::makeClusters( Extracter *e )
// void Merger::makeClusters( PathEntryMap *e )
// {
//     init(e);
//     makeSingletons();
// }

// void Merger::makeSingletons()
// {
//     PathEntryMap::const_iterator it;
//     for ( it = path_entries.begin(); it != path_entries.end(); ++it ) {
//         PathId pid = it->first;
//         AlignSummary aln;
//         aln.self(it->second.seq.size());
//         Cluster c;
//         c.add(pid, aln);
//         clusters.insert( std::pair<PathId, Cluster>( pid, c ) );
//     }
// }


void Merger::rebuildIndex()
{
    double t0 = mytime();
    if ( Param::verbose ) std::cout << "Rebuiding index with k-mer:" << kmer_size << "\n";
    kmer_paths.clear();
    for ( Clusters::iterator it = clusters.begin(); it != clusters.end(); ++it ) {
        PathId pid = it->first;
        if ( merged_paths.find(pid) != merged_paths.end() ) continue;
        std::string str = getPathSequence(pid);
        KmerArray  kmers = biostr::getKmers(str, kmer_size);
        kmer_paths.add( kmers, pid );
    }
    if ( Param::verbose >= 1 ) 
        std::cout << "Index rebuilt:" << mytime()-t0 << " sec\n";

    path_kmers.clear();
    path_kposs.clear();
    build();

    log.et_index += (mytime()-t0);
}


void Merger::makePathEntries( PathEntryMap & paths)
{
    PathEntryMap::iterator pit;

    for ( auto item : clusters ) {
        PathId pid = item.first;                                        
        if ( merged_paths.find(pid) != merged_paths.end() ) continue;
        
        //assert( id_map.find(pid) != id_map.end() );
        //PathId oid = id_map[pid];

        // pit = path_entries->find(oid);
        // assert( pit != path_entries->end() );

        //----------------
        // Update sequence
        //----------------
        Cluster c = item.second;
                    
        PathEntry e( c.sequence, c.lstop, c.rstop, c.ltrim, c.rtrim );
        paths.insert( std::pair<PathId, PathEntry>( pid, e ) );
    }
}

std::string Merger::getPathSequence(PathId &pid)
{
    // PathEntryMap::const_iterator it = path_entries.find(pid);
    // assert( it != path_entries.end() );
    
    // return it->second.seq;
    Clusters::const_iterator it = clusters.find(pid);
    assert( it != clusters.end() );
    return it->second.sequence;
}

//====================================================================
// Get the maximum numerical path ID from entries
//====================================================================
PathId Merger::getMaxPathId()
{
    PathId max = NOT_PATH;
    PathEntryMap::const_iterator it;
    for ( it = path_entries->begin(); it != path_entries->end(); ++it ) {
        PathId pid = it->first;
        assert( pid >= 0 );
        if ( pid > max ) max = pid;
    }
    return max;
}
