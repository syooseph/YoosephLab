#include "rjoiner.h"

ReadJoiner::ReadJoiner()
{

}

ReadJoiner::~ReadJoiner()
{

}


// BASE class (Connecter) pure function
void ReadJoiner::initMaps()
{

}

// BASE class (Connecter) pure function
void ReadJoiner::updateIndex( PathId spid,
                              std::string &old_pivot,
                              std::string &new_pivot,
                              int direction)
{

}

// BASE class (Connecter) pure function
bool ReadJoiner::latchOverlapPaths( PathId &spid, 
                                    int direction )
{
    return true;
}

// Function overloading
void ReadJoiner::run()
{

}


// void ReadJoiner::getGsaTypes( GsaTypeList &gsas,
//                               std::string &end_str )
// {
//     for ( int i = 0; i < Param::nparts; i++ ) {
//         int sr = i*int(nreads/Param::nparts);                
//         BoundType srch = gsa[i].search((SfaChar*)end_str.c_str(), end_str.size());
//         for ( int j = srch.first; j <= srch.second; j++ ) {
//             GsaType sa = gsa[i].getAt(j);
//             ReadId rid = sr + sa.doc;
//             if ( preads[rid] != NOT_PATH ) continue;
//             if ( preads[rid] == BAD_READ ) continue;

//             gsas.push_back( sa );
//         }
//     }
// }


// //====================================================================
// // Extract bridging reads between two paths
// //====================================================================
// bool ReadJoiner::extractSharedReads( const GsaTypeList *pivot_gsas,
//                                      const GsaTypeList *match_gsas,
//                                      std::tr1::unordered_map<ReadId,GsaType> &pivot_rmap, 
//                                      std::tr1::unordered_map<ReadId,GsaType> &match_rmap, 
//                                      std::tr1::unordered_map<ReadId, bool> &comm_reads,
//                                      PathId spid,
//                                      PathId mpid,
//                                      int direction )
// {
//     GsaTypeList::const_iterator gt;

//     //------------------------------------------
//     // First, save reads covering pivot sequence
//     //------------------------------------------
//     double tic = mytime();
//     for ( gt = pivot_gsas->begin(); gt != pivot_gsas->end(); ++gt ) {
//         if ( recruited.find(gt->doc) != recruited.end() ) continue;
//         pivot_rmap.insert( std::pair<ReadId,GsaType>( gt->doc, *gt ) );
//     }
//     bridge_log.et_read_pivot += (mytime()-tic);

//     if ( (int)pivot_rmap.size() < Param::min_share ) return false;

//     //------------------------------------------------
//     // Find reads covering match sequence.
//     // At the same time, find common reads from pivot.
//     //------------------------------------------------
//     tic = mytime();
//     for ( gt = match_gsas->begin(); gt != match_gsas->end(); ++gt ) {
//         if ( recruited.find(gt->doc) != recruited.end() ) continue;
//         if ( pivot_rmap.find(gt->doc) != pivot_rmap.end() ) {
//             comm_reads.insert( std::pair<ReadId,bool>(gt->doc, true) );
//             match_rmap.insert( std::pair<ReadId,GsaType>( gt->doc, *gt ) );
//         }
//     }
//     bridge_log.et_read_match += (mytime()-tic);

//     //----------------------------------------
//     // We need minimum number of common reads.
//     //----------------------------------------
//     if ( (int)comm_reads.size() < Param::min_share ) return false;

    
//     //--------------------------------------------------------------
//     // Extract reads that have valid aligned ranges againt two paths
//     //--------------------------------------------------------------
//     tic = mytime();
//     //ReadIdSet bad_reads;
//     std::unordered_map<ReadId, bool> bad_reads;
//     std::tr1::unordered_map<ReadId,GsaType>::iterator it;
//     for ( it = pivot_rmap.begin(); it != pivot_rmap.end(); ++it) {
//         ReadId rid = it->first;
//         if ( comm_reads.find(it->first) == comm_reads.end() ) {
//             //bad_reads.insert( rid );
//             bad_reads.insert( std::pair<PathId, bool>( rid, true) );
//             continue;
//         }
        
//         int pivot_rpos = it->second.pos;
//         int match_rpos = match_rmap[it->first].pos;

//         //---------------------------------------------------------
//         // Valid
//         //                   ooo..........  pivot
//         //           ooo.....ooo            read
//         // ..........ooo                    match
//         //           l       r
//         // Invalid
//         // We handled this already in long/short overlap extension
//         //           ooo.................  pivot
//         //           oooo.......           read
//         // ...........ooo                  match
//         //           rl
//         //--------------------------------------------------------
//         if ( direction == LEFT ) {
//             if ( match_rpos >= pivot_rpos ) {
//                 bad_reads.insert( std::pair<PathId, bool>( rid, true) );
//                 //bad_reads.insert( rid );
//                 continue;
//             } 
//         } 
//         //---------------------------------------------------------
//         // Valid
//         // ..........ooo                    pivot
//         //           ooo.....ooo            read
//         //                   ooo..........  match
//         //           l       r
//         // Invalid
//         // We handled this already in long/short overlap extension
//         // ...........ooo                  pivot
//         //           oooo.......           read
//         //           ooo.................  match
//         //           rl
//         //--------------------------------------------------------
//         else {
//             if ( match_rpos <= pivot_rpos ) {
//                 //bad_reads.insert( rid );
//                 bad_reads.insert( std::pair<PathId, bool>( rid, true) );
//                 continue;
//             }
//         }
//     }
//     bridge_log.et_read_bad += (mytime()-tic);    

//     tic = mytime();
//     for ( auto item : bad_reads ) {
//         pivot_rmap.erase(item.first);
//         comm_reads.erase(item.first);
//     }
//     bridge_log.et_read_trim += (mytime()-tic);

//     if ( (int)comm_reads.size() < Param::min_share ) return false;

//     return true;
// }

// //====================================================================
// // Drop reads with parital sequence match to pivot/match
// //====================================================================
// bool ReadJoiner::dropPartialReads( std::tr1::unordered_map<ReadId,GsaType> &pivot_rmap, 
//                                    std::tr1::unordered_map<ReadId,GsaType> &match_rmap, 
//                                    std::tr1::unordered_map<ReadId, bool> &comm_reads,
//                                    ReadPathPoss &pivot_poss, 
//                                    ReadPathPoss &match_poss,
//                                    std::string &pivot_str,
//                                    std::string &match_str,
//                                    int direction )
// {
//     // exact match of read
//     for ( auto it = pivot_rmap.begin(); it != pivot_rmap.end(); ) {
//         IntPair poss;
//         bool found = ( direction == LEFT ) ?
//             goodRead( poss, it->second, pivot_str, LEFT ) :
//             goodRead( poss, it->second, pivot_str, RIGHT ) ;
        
//         if ( found ) {
//             pivot_poss.insert( std::pair<ReadId, IntPair>( it->first, poss ) );
//             ++it;
//         } else {
//             comm_reads.erase(it->first);
//             match_rmap.erase(it->first);
//             pivot_rmap.erase(it++);
//         }
//     }

//     if ( (int)comm_reads.size() < Param::min_share ) return false;

//     // exact match of read
//     for ( auto it = match_rmap.begin(); it != match_rmap.end(); ) {
//         IntPair poss;
//         bool found = ( direction == LEFT ) ?
//             goodRead( poss, it->second, match_str, RIGHT ) :
//             goodRead( poss, it->second, match_str, LEFT ) ;

//         if ( found ) {
//             match_poss.insert( std::pair<ReadId, IntPair>( it->first, poss ) );
//             ++it;
//         } else {
//             comm_reads.erase(it->first);
//             pivot_poss.erase(it->first);
//             pivot_rmap.erase(it->first);
//             match_rmap.erase(it++);
//         }
//     }
    
//     if ( (int)comm_reads.size() < Param::min_share ) return false;

//     return true;
// }

// //====================================================================
// // Check prefix/suffix of read has identical sequence with pivot/match
// // We only know the read are in perfect match for the given length.
// // So, check the rest of sequence has perferct match or not.
// //==================================================================== 
// bool ReadJoiner::goodRead( IntPair &poss, 
//                           GsaType &sa, 
//                           std::string &pstr, 
//                           int direction )
// {
//     int l = Param::bridge_overlap;

//     ReadId r = sa.doc;
//     int    p = sa.pos;

//     assert( (int)r < nreads );
//     std::string rstr = seqs[r];

//     //-----------------------------
//     //           | l |
//     // read .....ooooo??
//     //           ooooo??...... path
//     //           p
//     //-----------------------------
//     if ( direction == LEFT ) {
//         std::string rsub = rstr.substr( p+l );
//         if ( rsub.size() > pstr.size()-l ) return false;
//         std::string psub = pstr.substr(l,rsub.size());
//         if ( rsub != psub ) return false;

//         // (read, path)
//         poss = IntPair( p, 0 );
//         return true;
//     } 
    
//     //-----------------------------
//     //             | l |
//     // path .....??ooooo
//     //           ??ooooo...... read
//     //             p
//     //-----------------------------
//     else {
//         std::string rsub = rstr.substr(0,p);
//         if ( rsub.size() > pstr.size()-l ) return false;
//         std::string psub = pstr.substr( pstr.size()-(l+rsub.size()), rsub.size() );
//         if ( rsub != psub ) return false;
        
//         // (read,path)
//         poss = IntPair(0, pstr.size()-(l+rsub.size()) );
        
//         return true;
//     }
        
// }

// //====================================================================
// // Now, find an evidence to stitch two paths.
// // Find a middle string which is either a string to connect two paths 
// // or tiny tiny overlap of two paths
// //====================================================================
// bool ReadJoiner::determineBridgingEvidence( ReadIdSet &bridges,
//                                            int &overlap,
//                                            std::string &mid_str,
//                                             std::string &pivot_str,
//                                             std::string &match_str,
//                                             ReadPathPoss &pivot_poss, 
//                                             ReadPathPoss &match_poss,
//                                             int direction )
// {
//     std::unordered_map< std::string, ReadIdSet > ladder_supports;
//     std::unordered_map< std::string, int > str_overlaps;

//     if ( Param::verbose ) {
//         std::cout << "pivot:" << pivot_str << "\n"
//         << "match:" << match_str << "\n"
//         << "direction:" << direction << "\n";
//     }

//     // get consistent middle string
//     for ( auto item : pivot_poss ) {
//         ReadId rid = item.first;
//         int pivot_rpos = item.second.first;
//         int pivot_ppos = item.second.second;
        
//         assert( match_poss.find(rid) != match_poss.end() );
//         int match_rpos = match_poss[rid].first;
//         int match_ppos = match_poss[rid].second;
        
//         std::string read_str = seqs[rid];

        
//         if ( Param::verbose ) std::cout << "rid:" << rid << "\t" << read_str << "\n";

//         int o = -1; // no overlap
//         //-------------------------------------------
//         //                     ooooo..........  pivot
//         //           ooooo.....ooooo            read
//         // ..........ooooo                      match
//         //           | l |     r
//         //-------------------------------------------
//         if ( direction == LEFT ) {
//             assert( match_ppos < (int)match_str.size() );
//             int l = (int)match_str.size()-match_ppos;
//             int r = pivot_rpos;

//             if ( Param::verbose ) 
//                 printf("l:%d, r:%d, direction:%d\n", l,r,direction);

//             // no overlap
//             if ( r >= l ) {
//                 mid_str = read_str.substr(l, r-l );
//             } else {
//                 o = abs(l-r);
//                 if ( Param::verbose ) std::cout << "overlap:" << o << "\n";
//                 assert( o < Param::bridge_overlap );
//                 mid_str = pivot_str.substr(0,o);
//             }
//             if ( Param::verbose ) 
//                 std::cout << "mid-str:" << mid_str << "\n";
//         } 
        
//         //-------------------------------------------
//         // ..........ooooo                      pivot
//         //           ooooo.....ooooo            read
//         //                     ooooo..........  match
//         //           | l |     r
//         //-------------------------------------------
//         else {
//             assert( pivot_ppos < (int)pivot_str.size() );
//             int l = (int)pivot_str.size()-pivot_ppos;
//             int r = match_rpos;

//             if ( Param::verbose ) 
//                 printf("l:%d, r:%d, direction:%d\n", l,r,direction);

//             // no overlap
//             if ( r >= l ) {
//                 mid_str = read_str.substr(l, r-l );
//             } else {
//                 o = abs(l-r);
//                 if ( Param::verbose ) std::cout << "overlap:" << o << "\n";
//                 assert( o < Param::bridge_overlap );
//                 mid_str = match_str.substr(0,o);

//             }

//             if ( Param::verbose ) 
//                 std::cout << "mid-str:" << mid_str << "\n";
//         } 

//         if ( ladder_supports.find( mid_str ) == ladder_supports.end() )
//             ladder_supports.insert( std::pair<std::string, ReadIdSet>( mid_str, ReadIdSet() ) );
//         ladder_supports[mid_str].insert(rid);

//         str_overlaps[mid_str] = o;
//     }
    
//     ReadIdSet   max_set;
//     std::string max_mid;
//     size_t max_num = 0;
//     for ( auto item : ladder_supports ) {
//         if ( item.second.size() > max_num ) {
//             max_num = item.second.size();
//             max_set = item.second;
//             max_mid = item.first;
//         }
//     }
    
//     if ( (int)max_num < Param::min_share ) return false;

//     mid_str = max_mid;
//     bridges = max_set;
//     overlap = str_overlaps[max_mid];
//     return true;    
// }




void ReadJoiner::updateSequence( PathId spid,
                                 PathId mpid,
                                 int overlap,
                                 std::string &mid_str,
                                 int direction )
{
    assert( merged_paths.find(spid) == merged_paths.end() );
    assert( merged_paths.find(mpid) == merged_paths.end() );

    std::string sstr = extenders[spid].sequence;
    std::string mstr = extenders[mpid].sequence;
    
    std::string nstr;
    if ( direction == LEFT ) {
        if ( overlap > 0 ) 
            nstr = mstr + sstr.substr(overlap);
        else
            nstr = mstr + mid_str + sstr;
    } else {
        if ( overlap > 0 ) 
            nstr = sstr + mstr.substr(overlap);
        else
            nstr = sstr + mid_str + mstr;
    }
    extenders[spid].sequence = nstr;

    if ( Param::verbose ) std::cout << "nstr:" << nstr << "\n";
}

void ReadJoiner::adjustAlignPositions( AlignSummary &summary, 
                                       const int &overlap, 
                                       const std::string &mid_str, 
                                       const PathId &spid, 
                                       const PathId &mpid, 
                                       const int &direction )
{
    if ( direction == LEFT ) {
        summary.lgap.first = ( overlap > 0 ) ?
            extenders[mpid].sequence.size()-overlap :
            extenders[mpid].sequence.size()+mid_str.size();
        summary.lgap.second = 0;
    } else {
        summary.lgap.first = 0;
        summary.lgap.second = ( overlap > 0 ) ?
            extenders[spid].sequence.size()-overlap :
            extenders[spid].sequence.size()+mid_str.size();
    }
}

// void ReadJoiner::addReads( PathId spid,
//                            PathId mpid,
//                            ReadIdSet &bridges,
//                            ReadPathPoss &pivot_poss,
//                            ReadPathPoss &match_poss,
//                            AlignSummary &summary,
//                            int direction )
// {
//     if ( added_reads.find(mpid) != added_reads.end() ) {
//         if ( Param::verbose ) std::cout << "Existing reads for mpid:" << mpid << "\tsize:" << added_reads[mpid].size() << "\n";
//         ReadEntryList::iterator jt;
//         for ( jt = added_reads[mpid].begin(); jt != added_reads[mpid].end(); ++jt ) {
//             if ( direction == RIGHT ) 
//                 jt->ppos += summary.lgap.second;
//             added_reads[spid].push_back(*jt);
//         }
//         added_reads.erase(mpid);
//     }

//     if ( Param::verbose ) std::cout << "New reads for spid:" << spid << "\tsize:" << bridges.size() << "\n";
//     for ( auto rid : bridges ) {
        
//         unsigned rpos, ppos;
//         if ( direction == LEFT ) {
//             assert( match_poss.find(rid) != match_poss.end() );
//             rpos = match_poss[rid].first;
//             ppos = match_poss[rid].second;
//         } else {
//             assert( pivot_poss.find(rid) != pivot_poss.end() );
//             rpos = pivot_poss[rid].first;
//             ppos = pivot_poss[rid].second;
//         }

//         if ( Param::verbose ) std::cout << "rid:" << rid << "\trpos:" << rpos << "\tppos:" << ppos << "\n";
//         added_reads[spid].push_back( ReadEntry(rid,rpos,ppos) );

//         recruited.insert( std::pair<ReadId, PathId>(rid,spid) );
//     }

// }


void ReadJoiner::printAddedReads()
{
    std::cout << "Added reads\t#Paths:" << added_reads.size() << "\n";
    PathReadsMap::iterator jt;
    size_t i;
    for ( i = 1, jt = added_reads.begin(); jt != added_reads.end(); ++i, ++jt ) {
        std::cout << i << "\tPathId:" << jt->first << "\t";
        std::cout << "# reads:" << jt->second.size() << "\n";
        ReadEntryList::iterator rt;
        for ( rt = jt->second.begin(); rt != jt->second.end(); ++rt )
            std::cout << rt->read << "\t" << rt->rpos << "\t" << rt->ppos << "\n";
    }
}

size_t ReadJoiner::getAddedReadsCount()
{
    size_t count = 0;
    for ( auto item : added_reads )
        count += item.second.size();
    return count;
}


void ReadJoiner::purgeTemp()
{

}
