#include "ljoiner.h"

LongJoiner::LongJoiner()
{

}

LongJoiner::~LongJoiner()
{

}

void LongJoiner::initMaps()
{
    
    for ( ExtenderMap::iterator it = extenders.begin(); it != extenders.end(); ++it ) {
        PathId pid = it->first;
        std::string str = it->second.sequence;
        if ( (int)str.size() < Param::extend_length ) continue;
        addPathIds( pid, str );
    }
}

void LongJoiner::purgeTemp()
{
    path_map.clear();
    path_kmers.clear();
}


bool LongJoiner::latchableCandidate( int &qs, 
                                     int &ss, 
                                     int &qe, 
                                     int &se,
                                     PathId spid, 
                                     PathId mpid, 
                                     std::string &pivot, 
                                     std::string &match, 
                                     const KmerArray *pkmers,
                                     const KmerArray *mkmers,
                                     const KmerPosMap *pposs,
                                     const KmerPosMap *mposs,
                                     int direction )
{
    if ( direction == LEFT && match[match.size()-1] == alpha::STOP_CODON ) {
        if ( !Param::ignore_stop ) return false;
    }

    if ( Param::check_stop && stopConflict(spid, mpid, direction) ) {
        if ( Param::verbose ) std::cout << "Stop conflicts\n";
        return false;
    }

    //--------------------------------------------------------
    // Find a latchable sequence range between pivot and match
    //--------------------------------------------------------
    double tic = mytime();
    LatchSectioner sect( &match, &pivot, mkmers, pkmers, mposs, pposs );
    sect.setDirection( direction );    
    bool found = sect.find();
    log.et_section += (mytime()-tic);
    if ( Param::verbose ) std::cout << "section time:" << mytime()-tic << " sec\n";
    if ( !found ) return false;
    
    Section s = sect.getSection();
    ss = s.sbeg; se = s.send;
    qs = s.qbeg; qe = s.qend;
    
    return true;
}

// void LongJoiner::makePosPairMap( std::multimap<int,int> &poss_pair,
//                                 KmerArray &query_kmers,
//                                 KmerArray &sbjct_kmers )
// {
//     std::tr1::unordered_map<KmerId, std::vector<int> > qposs, sposs;
//     for ( size_t i = 0; i < query_kmers.size(); i++ ) 
//         qposs[query_kmers[i]].push_back(i);
//     for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) 
//         sposs[sbjct_kmers[i]].push_back(i);


//     std::tr1::unordered_map<KmerId, std::vector<int> >::iterator it;
//     //std::multimap<int,int> poss_pair;
//     for ( it = sposs.begin(); it != sposs.end(); ++it ) {
//         std::vector<int> sp = it->second;
//         if ( qposs.find(it->first) == qposs.end() ) continue;
//         std::vector<int> qp = qposs[it->first];
//         for ( size_t i = 0; i < sp.size(); i++ ) 
//             for ( size_t j = 0; j < qp.size(); j++ )
//                 poss_pair.insert(IntPair(qp[j], sp[i]));
//     }
// }

// bool LongJoiner::findRegion( std::multimap<int,int> &poss_pair,
//                             std::multimap<int,int>::iterator &pt,
//                             int &qs,
//                             int &ss,
//                             int &qe, 
//                             int &se,
//                             int &ct,
//                             KmerArray &query_kmers )
// {
//     qs = qe = pt->first;  // query
//     ss = se = pt->second; // sbjct
//     ct = 0;
    
//    // if ( Param::verbose ) {
//    //      std::cout << "Kmer:" << alpha::IntegerToAminoAcid(query_kmers[qs], filter_kmer) << "\n";
//    //      printf("init\tqs:%d ss:%d qe:%d se:%d\n", qs,ss,qe,se);
//    //  }

//     std::multimap<int,int>::iterator qt;
//     for ( qt = pt; qt != poss_pair.end(); ++qt ) {
//         // if ( Param::verbose ) {
//         //     std::cout << "Kmer:" << alpha::IntegerToAminoAcid(query_kmers[qt->first], filter_kmer) << "\n";
//         //     printf("loop\tqs:%d ss:%d qe:%d se:%d\n", qs,ss, qt->first, qt->second);
//         // }
//         if ( qt->first  < qs  || qt->first  < qe ) continue;//break;
//         if ( qt->second < ss  || qt->second < se ) continue;//break;
//         int diff = (qt->first-qs) - (qt->second-ss);
//         if ( abs(diff) > Param::extend_off_nbase ) continue;//break;
//         qe = qt->first;
//         se = qt->second;
//         ct++;
//     }
    
//     if ( Param::verbose ) printf("region:\tqs:%d ss:%d qe:%d se:%d\tmatch:%d\n", qs,ss, qe,se, ct);
//     int diff = (qe-qs) - (se-ss);
//     //if ( Param::verbose ) std::cout << "diff:" << diff << "\n";
//     //if ( abs(diff) > 1 ) return false;//continue;
//     if ( abs(diff) > Param::extend_off_nbase ) return false;//continue;
//     if ( ct < mink ) return false;//continue;
    
//     return true;
// }


// bool LongJoiner::extendRegion( int &qs,
//                               int &ss,
//                               int &qe,
//                               int &se,
//                               int qnkmers,
//                               int snkmers,
//                               int direction )
// {
//     //             ss--------------------se (sbjct)
//     // qs--------------------qe (query)
//     if ( direction == LEFT ) {
//         if ( ss > qs ) return false;
//         if ( ss > 0 ) {
//             qs -= ss; ss =0;
//         }
//         int off = qnkmers - (qe+1);
//         if ( off > 0 && se < snkmers-1 ) {
//             se += off;
//             if ( se >= snkmers ) return false;
//             qe = qnkmers - 1;
//         }
//     } 
//     // ss--------------------se (sbjct)
//     //             qs--------------------qe (query)
//     else {
//         if ( qs > ss ) return false;
//         if ( qs > 0 ) {
//             ss -= qs; qs =0;
//         }

//         int off = snkmers - (se+1);
//         if ( off > 0 && qe < qnkmers-1 ) {
//             qe += off;
//             if ( qe >= qnkmers ) return false;
//             se = snkmers -1;
//         }
//     }
//     return true;
//     //if ( Param::verbose ) printf("extended:\tqs:%d ss:%d qe:%d se:%d\n", qs,ss,qe,se);
// }



// bool LongJoiner::passKmerFilter( int qs,
//                                 int ss,
//                                 int qe, 
//                                 int se,
//                                 KmerArray &query_kmers,
//                                 KmerArray &sbjct_kmers )
// {
//     int len_ext = (qe-qs < se-ss) ? 
//         (qe-qs+filter_kmer) :
//         (se-ss+filter_kmer) ;    
//     assert(len_ext>0);

//     int min_ext = filter::minSameKmerCount( len_ext, filter_kmer, 1-filter_score ) ;

//     int num_ext = 0;
//     //int a,b;
    
//     std::set<KmerId> query_set;
//     for ( int a = qs; a <= qe; a++ ) 
//         query_set.insert( query_kmers[a] );
//     for ( int b = ss; b <= se; b++ )
//         if ( query_set.find(sbjct_kmers[b]) != query_set.end() )
//             num_ext++;
    
//     if ( Param::verbose ) printf("extended region:\tqs:%d ss:%d qe:%d se:%d\n", qs,ss,qe,se);
//     if ( Param::verbose ) std::cout << "extended length:" << len_ext << "\tmin k-mer match:" << min_ext << "\t#match:" << num_ext << "\n";

//     if ( num_ext < min_ext ) return false;
//     return true;        
// }

// make index of new sequence
void LongJoiner::addPathIds( //KmerToPathMap &path_map,
                            PathId pid,
                            std::string sequence )
{
    // kmers for new sequence
    KmerArray kmers = biostr::getKmers(sequence, Param::extend_anchor_kmer);
    for ( size_t i = 0; i < kmers.size(); i++ )
        path_map[kmers[i]].insert(pid);


    path_kmers.add(pid, kmers);

    KmerPosMap kposs;
    for ( size_t j = 0; j < kmers.size(); j++ )
        kposs[kmers[j]].push_back(j);
    path_kposs.add(pid, kposs);
}


// drop entries of old sequence
void LongJoiner::dropPathIds( //KmerToPathMap &path_map,
                             PathId pid,
                             std::string sequence )
{
    //KmerArray kmers = biostr::getKmers(sequence, Param::extend_anchor_kmer);
    const KmerArray *kmers = path_kmers.getPointer(pid);
    for ( size_t i = 0; i < kmers->size(); i++ )
        path_map[(*kmers)[i]].erase(pid);

    path_kmers.erase(pid);
    path_kposs.erase(pid);
}

bool LongJoiner::latchOverlapPaths(  PathId &spid, int direction )
{
    std::string pivot = extenders[spid].sequence;
    
    double tic = mytime();
    std::vector<JoinEntry> entries;
    findPathsByKmerCount( entries, spid, direction );
    log.et_path += (mytime()-tic);
    log.ct_path += entries.size();
    if ( entries.size() == 0 ) return false;

    //KmerArray pkmers = biostr::getKmers(pivot, Param::extend_anchor_kmer);
    const KmerArray  *pkmers = path_kmers.getPointer(spid);
    const KmerPosMap *pposs = path_kposs.getPointer(spid);
    // KmerPossMap pposs;      ///< pivot kmer positions
    // for ( size_t i = 0; i < pkmers->size(); i++ ) 
    //     pposs[(*pkmers)[i]].push_back(i);

    tic = mytime();
    /* Find reads in library size length */
    std::tr1::unordered_map<ReadId, bool> pivot_reads;
    if ( Param::pair_flag && Param::ljoin_pe_support ) {
        std::string sub_pivot = getStringInRange(pivot, direction, true );            
        if ( Param::verbose ) {
            std::cout << "sub-pivot:" << sub_pivot << "\n";
            std::cout << "length:" << sub_pivot.size() << "\n";
        }
        findReads( sub_pivot, pivot_reads );
    }
    log.et_reads += (mytime()-tic);
    if ( Param::verbose ) 
        printf("# pivot reads in library range:%zu (%.4f sec)\n", pivot_reads.size(), mytime()-tic );
    
    tic = mytime();    
    bool success = false;
    size_t i;
    for ( i = 0; i < entries.size(); i++ ) {
        //PathId mpid = entries[i].pid;
        log.ct_trial++;
        success = tryConnectOverlap(entries[i], spid, pivot, pkmers, pposs, pivot_reads, direction );
        if ( success ) {
            log.ct_latch_succ++;
            break;
        }
    }
    log.et_latch += (mytime()-tic);

    if ( Param::verbose ) 
        printf("#entries:%zu\titer:%zu\tsuccess:%d\ttime:%.4f\n", entries.size(), i, success, mytime()-tic);
                                 

    return success;
}


void LongJoiner::findPathsByKmerCount( std::vector<JoinEntry> &entries, 
                                        PathId spid, 
                                        int direction )
{
    double tic = mytime();
    std::map<PathId, int> count_map;
    getKmerCountMap( count_map, spid, direction );
    log.et_path_filter += (mytime()-tic);

    if ( Param::verbose ) std::cout << "# candidates:" << count_map.size() << "\t" << mytime()-tic << " sec\n";
    if ( count_map.size() == 0 ) return;

    tic = mytime();
    entries.reserve(count_map.size());
    for ( std::map<PathId, int>::iterator ct =  count_map.begin(); ct != count_map.end(); ++ct ) {
        PathId pid = ct->first;
        int num = ct->second;
        //if ( num < mink ) continue;
        if ( num < Param::extend_anchor_mink ) continue;
        std::string seq = extenders[pid].sequence;
        entries.push_back( JoinEntry(pid, num, seq.size() ) );
    }
    log.et_path_entry += (mytime()-tic);

    tic = mytime();
    std::sort( entries.begin(), entries.end(), cmp_entry );
    log.et_path_sort += (mytime()-tic);
    
    if ( Param::verbose ) std::cout << "# good candidates:" << entries.size() << "\n";
    if ( Param::verbose ) {
        for ( size_t i = 0; i < entries.size(); i++ )
            std::cout << "pid:" << entries[i].pid 
                      << "\tmatch-kmer:" << entries[i].nkmers 
                      << "\tsequence-length:" << entries[i].seqlen << "\n";
    }
    //if ( entries.size() == 0 ) return;
}


bool LongJoiner::tryConnectOverlap(JoinEntry &entry, 
                                   PathId &spid, 
                                   std::string &pivot, 
                                   const KmerArray *pkmers, 
                                   const KmerPosMap *pposs,
                                   std::tr1::unordered_map<ReadId, bool> &pivot_reads,
                                   int direction )
{
    PathId mpid = entry.pid;
    int numk = entry.nkmers;
    std::string match = extenders[mpid].sequence;
    assert( match.size() > 0 );

    if ( Param::verbose ) {
        std::cout << "Pid:" << mpid << "\t#Filter-kmers:" << numk << "\tSequence-length:" << match.size() << "\n";
        std::cout << "Match:" << match << "\n";
    }
    
    //KmerArray mkmers = biostr::getKmers(match, Param::extend_anchor_kmer);
    const KmerArray* mkmers = path_kmers.getPointer(mpid);
    if ( mkmers->size() == 0 ) return false;

    const KmerPosMap *mposs = path_kposs.getPointer(mpid);
    // KmerPossMap mposs;
    // for ( size_t i = 0; i < mkmers->size(); i++ ) 
    //     mposs[(*mkmers)[i]].push_back(i);

    double tic = mytime();
    int match_beg, pivot_beg, match_end, pivot_end;
    //Section s;
    bool latchable = latchableCandidate(match_beg, pivot_beg, match_end, pivot_end, spid, mpid, pivot, match, pkmers, mkmers, pposs, mposs, direction);
    log.et_latch_check += (mytime()-tic);
    if ( !latchable ) return false;
    
    if ( direction == RIGHT ) {
        if ( match_beg != 0 ) {
            if ( Param::verbose) std::cout << "Wrong match region start\n";
            return false;
        }
    } else {
        if ( pivot_beg != 0 ) {
            if ( Param::verbose ) std::cout << "Wrong pivot region start\n";
            return false;
        }
    }

    // std::string sub_pivot = pivot.substr( pivot_beg, pivot_end-pivot_beg+Param::kmer_size);
    // std::string sub_match = match.substr( match_beg, match_end-match_beg+Param::kmer_size);
    std::string sub_pivot = pivot.substr( pivot_beg, pivot_end-pivot_beg+1);
    std::string sub_match = match.substr( match_beg, match_end-match_beg+1);

    if ( Param::verbose ) {
        printf("Pivot sub:\ts:%d\te:%d\n", pivot_beg, pivot_end );
        std::cout << sub_pivot << "\n";
        printf("Match sub:\ts:%d\te:%d\n", match_beg, match_end );
        std::cout << sub_match << "\n";
    }


    AlignSummary summary;
    
    bool good = false;
    if ( Param::align_base_first ) 
        good = checkByBases( summary, sub_pivot, sub_match );

    if ( !good )
        good = checkByAlign( summary, sub_pivot, sub_match, direction );
    
    if ( !good ) return false;
    

    int nlen =  pivot.size()+match.size()-summary.length;
    if ( Param::verbose ) std::cout << "New sequence size:" << nlen << "\n";
    if ( Param::pair_flag && Param::ljoin_pe_support ) {
        if ( nlen >= min_libsize ) {
            double lt0 = mytime();
            size_t npairs = findPairendReads( pivot_reads, pivot, match, direction );
            log.et_reads += (mytime()-lt0);
            if ( Param::verbose ) 
                printf("# match pairs in library range:%zu (%.4f sec)\n", npairs, mytime()-lt0 );
            
            if ( (int)npairs < Param::min_pair_reads ) {
                if ( Param::verbose ) std::cout << "No read support\n";
                log.ct_fail_pread++;
                return false; //continue;
            }
        }
        else {
            if ( Param::verbose ) 
                std::cout << "Too short length for pair-end support check\n";
            log.ct_fail_short++;
            return false;
        }
    }

    if ( direction == RIGHT && pivot_beg > 0 ) {
        summary.lgap.second += pivot_beg;
        summary.shift(pivot_beg, true);
        // if ( match_end+Param::kmer_size < (int)match.size() ) 
        //     summary.egap.first += ( match.size()- (match_end+Param::kmer_size) );
        if ( match_end+1 < (int)match.size() ) 
            summary.egap.first += ( match.size()- (match_end+1) );
    }
    if ( direction == LEFT && match_beg > 0 ) {
        summary.lgap.first += match_beg;
        summary.shift(match_beg, false);
        // if ( pivot_end+Param::kmer_size < (int)pivot.size() ) 
        //     summary.egap.second += ( pivot.size()- (pivot_end+Param::kmer_size) );
        if ( pivot_end+1 < (int)pivot.size() ) 
            summary.egap.second += ( pivot.size()- (pivot_end+1) );
    }

    if ( Param::verbose ) summary.print(std::cout);

    connect(spid, mpid, pivot, match, summary, direction);

    return true;
}

bool LongJoiner::checkByBases( AlignSummary &summary,
                               std::string &sub_pivot,
                               std::string &sub_match )
{
    if ( sub_pivot.size() != sub_match.size() ) 
        return false;

    double t0 = mytime();
    compareBases( summary, sub_pivot, sub_match );
    log.et_bases += (mytime()-t0);
    log.ct_bases ++;
    
    if ( summary.length < Param::extend_length ) 
        return false;

    if ( summary.posrate < Param::extend_score ) 
        return false;

    log.ct_bases_succ++;
    return true;
}

bool LongJoiner::checkByAlign( AlignSummary &summary,
                               std::string &sub_pivot,
                               std::string &sub_match,
                               int direction )
{
    double tic = mytime();

    summary.init();
    alignPair(summary, sub_pivot, sub_match, direction);
    log.et_align += (mytime()-tic);
    log.ct_align ++;

    if ( summary.lgap.first > 0 && summary.lgap.second > 0 ) return false;
    if ( summary.egap.first > 0 && summary.egap.second > 0 ) return false;
    if ( direction == LEFT  && summary.lgap.second > 0 ) return false;
    if ( direction == RIGHT && summary.lgap.first  > 0 ) return false; //continue;
    
    double score = !Param::identity_flag ? 
        summary.posrate :
        (double)summary.match/summary.length;


    if ( Param::verbose ) std::cout << "score:" << score << "\n";
    if ( score < Param::extend_score ) return false; //continue;
    if ( summary.length < Param::extend_length ) {
        if ( Param::verbose ) std::cout << "Short overlap\n"; 
        return false;
    }
    
    log.ct_align_succ++;    
    if (summary.ilist.size() || summary.dlist.size() ) log.ct_align_succ_indel++;

    return true;

}

void LongJoiner::updateIndex( PathId spid,
                             std::string &old_pivot,
                             std::string &new_pivot,
                             int direction )
{
    dropPathIds( spid, old_pivot );
    addPathIds ( spid, new_pivot );
}

//std::map<PathId, int>  
void LongJoiner::getKmerCountMap( std::map<PathId, int> &count_map,
                                 //PathIdSet skip_pids,
                                 PathId spid,
                                 int direction )
{
    //std::map<PathId, int> count_map;
    std::string pivot = extenders[spid].sequence;
    if ( (int)pivot.size() < Param::extend_length ) return;//return count_map;

    // std::string sub_pivot = direction==LEFT ?
    //     pivot.substr(0, pivot.size()/2) :
    //     pivot.substr(pivot.size()/2);
    std::string sub_pivot = direction==LEFT ?
        pivot.substr(0, Param::extend_length) :
        pivot.substr(pivot.size()-Param::extend_length);
    
    if ( Param::verbose) std::cout << "Pivot end-str:" << sub_pivot << "\n";
    //KmerArray kmers = biostr::getKmers(sub_pivot, Param::extend_anchor_kmer);

    const KmerArray *all_kmers = path_kmers.getPointer(spid);    
    int ncheck = Param::extend_length - Param::extend_anchor_kmer + 1;
    KmerArray kmers = direction == LEFT ? 
        KmerArray( all_kmers->begin(),  all_kmers->begin()  + ncheck ) :
        KmerArray( all_kmers->rbegin(), all_kmers->rbegin() + ncheck) ;
    
    KmerSet kset = KmerSet( kmers.begin(), kmers.end() );

    // if ( Param::verbose) {
    //     for ( KmerSet::iterator it = kset.begin(); it != kset.end(); ++it )
    //         std::cout << alpha::IntegerToAminoAcid(*it, Param::extend_anchor_kmer) << " ";
    //     std::cout << "\n";
    // }

    //KmerArray kmers = biostr::getKmers(pivot, Param::extend_filter_kmer);
    //KmerSet kset = KmerSet( kmers.begin(), kmers.end() );
    for ( KmerSet::iterator kit = kset.begin(); kit != kset.end(); ++kit ) {
        std::set<PathId>::iterator pit;
        //if ( !short_flag ) {
            if ( path_map.find(*kit) == path_map.end() ) continue;
            for ( pit = path_map[*kit].begin(); pit != path_map[*kit].end(); ++pit ) {
                if ( *pit == spid ) continue;
                if ( merged_paths.find(*pit) != merged_paths.end() ) continue; 
                if ( count_map.find(*pit) == count_map.end() )
                    count_map.insert( std::pair<PathId, int>(*pit, 0) );
                count_map[*pit]++;
            }
            //} 
        // else {
        //     if ( direction == LEFT ) {
        //         if ( rpath_map.find(*kit) == rpath_map.end() ) continue;
        //         for ( pit = rpath_map[*kit].begin(); pit != rpath_map[*kit].end(); ++pit ) {
        //             //if ( skip_pids.find(*pit) != skip_pids.end() ) continue; 
        //             if ( *pit == spid ) continue;
        //             if ( merged_paths.find(*pit) != merged_paths.end() ) continue; 
        //             if ( count_map.find(*pit) == count_map.end() )
        //                 count_map.insert( std::pair<PathId, int>(*pit, 0) );
        //             count_map[*pit]++;
        //         }
        //     } else {
        //         if ( lpath_map.find(*kit) == lpath_map.end() ) continue;
        //         for ( pit = lpath_map[*kit].begin(); pit != lpath_map[*kit].end(); ++pit ) {
        //             //if ( skip_pids.find(*pit) != skip_pids.end() ) continue; 
        //             if ( *pit == spid ) continue;
        //             if ( merged_paths.find(*pit) != merged_paths.end() ) continue; 
        //             if ( count_map.find(*pit) == count_map.end() )
        //                 count_map.insert( std::pair<PathId, int>(*pit, 0) );
        //             count_map[*pit]++;
        //         }
        //     }
        // }
    }
    //return count_map;
}
