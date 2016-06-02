#include "sjoiner.h"

TinyJoiner::TinyJoiner()
{

}

TinyJoiner::~TinyJoiner()
{

}

void TinyJoiner::initMaps()
{
    
    for ( ExtenderMap::iterator it = extenders.begin(); it != extenders.end(); ++it ) {
        PathId pid = it->first;
        std::string str = it->second.sequence;
        if ( (int)str.size() < Param::extend_length ) continue;
        updateEndString(str, pid);                    
    }
}

void TinyJoiner::purgeTemp()
{
    lend_map.clear();
    rend_map.clear();
}

void TinyJoiner::updateEndString( std::string &str,
                                 PathId pid )
{
    for ( int i = Param::suffix_overlap; i < Param::extend_length; i++ ) {
        if ( i >= (int)str.size() ) break;
        std::string sub = str.substr(0, i);
        lend_map[i][sub].insert(pid);
    }

    for ( int i = Param::suffix_overlap; i < Param::extend_length; i++ ) {
        if ( i >= (int)str.size() ) break;
        std::string sub = str.substr(str.size()-i);
        rend_map[i][sub].insert(pid);
    }
}

void TinyJoiner::deleteEndString( std::string &str,
                                 PathId pid )
{
    for ( int i = Param::suffix_overlap; i < Param::extend_length; i++ ) {
        if ( i >= (int)str.size() ) break;
        std::string sub = str.substr(0, i);
        lend_map[i][sub].erase(pid);
    }

    for ( int i = Param::suffix_overlap; i < Param::extend_length; i++ ) {
        if ( i >= (int)str.size() ) break;
        std::string sub = str.substr(str.size()-i);
        rend_map[i][sub].erase(pid);
    }
}



bool TinyJoiner::latchOverlapPaths( PathId &spid, int direction )
{
    std::tr1::unordered_map<PathId, bool> checked;

    std::string pivot = extenders[spid].sequence;

    double tic = mytime();
    std::vector<size_t> lengths;
    std::vector<PathId> matches;
    findPathsByEndString( matches, lengths, spid, pivot, checked, direction ) ;
    log.et_path += (mytime()-tic);

    if ( lengths.size() == 0 ) return false;
    
    tic = mytime();
    /* Find reads in library size length */
    std::tr1::unordered_map<ReadId, bool> pivot_reads;
    if ( Param::pair_flag && Param::sjoin_pe_support ) {
        std::string sub_pivot = getStringInRange(pivot, direction, true );            
        if ( Param::verbose ) std::cout << "Sub-pivot in range:" << sub_pivot << "\tlength:" << sub_pivot.size() << "\n";
        findReads( sub_pivot, pivot_reads );
    }
    log.et_reads += (mytime()-tic);
    if ( Param::verbose ) 
        printf("# pivot reads in library range:%zu (%.4f sec)\n", pivot_reads.size(), mytime()-tic );

    tic = mytime();
    bool success = false;
    size_t i;
    for ( i = 0; i < matches.size(); i++ ) {
        if ( Param::check_stop && stopConflict(spid, matches[i], direction) ) {
            if ( Param::verbose ) std::cout << "Stop conflicts\n";
            continue;
        }
        
        log.ct_trial++;
        success = tryConnectOverlap( lengths[i], matches[i], spid, pivot, pivot_reads, direction );
        if ( success ) {
            log.ct_latch_succ++; 
            //return true;
            break;
        }
        else checked.insert(std::pair<PathId,bool>(matches[i], true));
    }
    log.et_latch += (mytime()-tic);

    if ( Param::verbose ) 
        printf("#matches:%zu\titer:%zu\tsuccess:%d\ttime:%.4f\n", matches.size(), i, success, mytime()-tic);

    return success;
}

bool TinyJoiner::tryConnectOverlap( size_t length,
                                    PathId &mpid,
                                    PathId &spid, 
                                    std::string &pivot, 
                                    std::tr1::unordered_map<ReadId, bool> &pivot_reads,
                                    int direction )
{
    std::string opivot = pivot;

    std::string match = extenders[mpid].sequence;
    if ( Param::verbose ) {
        std::cout << "Pivot:" << pivot << "\n";
        std::cout << "Direction:" << direction << "\n";
        std::cout << "Match Pid:" << mpid << "\t#length:" << length << "\tSequence-length:" << match.size() << "\n";
        std::cout << "Match:" << match << "\n";
    }
    
    std::string ladder;
    if ( ! getLadderString( ladder, length, pivot, match, direction ) ) return false;
    
    if ( Param::verbose ) std::cout << "Ladder:" << ladder << "\n";
    
    int sum = 0;
    for ( int i = 0; i < Param::nparts; i++ ) {
        BoundType srch = gsa[i].search((SfaChar*)ladder.c_str(), ladder.size());
        if ( Param::verbose ) 
            std::cout << "(" << srch.first << "\t" << srch.second << ")\n";
        
        int count = srch.second>=srch.first ? (srch.second-srch.first+1) : 0;
        sum += count;
    }

    if ( Param::verbose ) {
        std::cout << "#SFA search:" << sum << "\n";
        if ( sum == 0 ) std::cout << "Invalid ladder string\n";
        //else std::cout << "Short overlap extension success\n";
    }
    if ( sum == 0 ) return false;
    
    int nlen =  pivot.size()+match.size()-length;
    //if ( Param::verbose ) std::cout << "New sequence size:" << nlen << "\n";
    if ( Param::pair_flag && Param::sjoin_pe_support ) {
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
    
    AlignSummary summary;
    summary.score = 1;
    summary.length = summary.match = length;
    summary.posrate = 1.0;
    if ( direction == LEFT ) {
        summary.lgap.first  = match.size()-length;
        summary.lgap.second = 0;
        summary.egap.first  = 0;
        summary.egap.second = pivot.size()-length;
    } else {
        summary.egap.first  = match.size()-length;
        summary.egap.second = 0;
        summary.lgap.first  = 0;
        summary.lgap.second = pivot.size()-length;
    }

    if ( Param::verbose ) summary.print(std::cout);

    connect(spid, mpid, pivot, match, summary, direction);

    //std::string npivot = extenders[spid].sequence;
    //deleteEndString(opivot, spid);
    //updateEndString(npivot, spid);        
    
    return true;
}


bool TinyJoiner::getLadderString( std::string &ladder,
                                 size_t length,
                                 std::string &pivot,
                                 std::string &match,
                                 int direction )
{
    //std::string lstr, rstr;
    std::string sub_pivot, sub_match;
    if ( direction == LEFT ) {
        sub_pivot = pivot.substr(0,length);
        sub_match = match.substr(match.size()-length);
    } else {
        sub_pivot = pivot.substr(pivot.size()-length);
        sub_match = match.substr(0,length);
    }

    if ( Param::verbose ) {
        std::cout << "sub-pivot:" << sub_pivot << "\n";
        std::cout << "sub-match:" << sub_match << "\n";
    }

    assert(sub_pivot==sub_match);


    int half = (Param::back_trace - (int)length ) / 2;
    if ( (int)match.size() < half+(int)length ||
         (int)pivot.size() < half+(int)length ) return false;
    
    if ( direction == LEFT ) {
        ladder = match.substr( match.size()-(half+length) );
        ladder += pivot.substr( length, half );
    } else {
        ladder = pivot.substr( pivot.size()-(half+length) );
        ladder += match.substr( length, half );
    }
    return true;
}



/** 
 * Only shorter sequence is allowed because of short latch region
 */
void TinyJoiner::findPathsByEndString( std::vector<PathId> &matches,
                                      std::vector<size_t> &lengths, 
                                      PathId spid, 
                                      std::string &pivot,
                                      std::tr1::unordered_map<PathId, bool> checked,
                                      int direction )
{
    EndStringMap *end_map = direction == LEFT ? &rend_map : &lend_map;
    
    for ( int i = Param::extend_length; i>=Param::suffix_overlap; i-- ) {
        if ( (int)pivot.size() <= i ) continue;

        std::string tip = direction == LEFT ? 
            pivot.substr(0,i) : pivot.substr(pivot.size()-i);
        if ( Param::verbose ) std::cout << "Tip-str:" << tip << "\n";
        if ( end_map->find(i) == end_map->end() ) continue;
        if ( (*end_map)[i].find(tip) == (*end_map)[i].end() )  continue;

        std::set<PathId>::iterator it;
        for ( it = (*end_map)[i][tip].begin(); it != (*end_map)[i][tip].end(); ++it ) {
            if ( *it == spid ) continue;
            if ( merged_paths.find(*it) != merged_paths.end() ) continue;
            if ( checked.find(*it) != checked.end() ) continue;

            // Length contraint //
            std::string match = extenders[*it].sequence;
            if ( match.size() > pivot.size() ) continue;
            
            if ( Param::verbose ) 
                std::cout << "Found:\tpath:" << *it << "\tlength:" << i << "\n";
            matches.push_back(*it);
            lengths.push_back(i);
        }
    }
}


void TinyJoiner::updateIndex( PathId spid,
                             std::string &old_pivot,
                             std::string &new_pivot,
                             int direction )
{
    deleteEndString(old_pivot, spid);
    updateEndString(new_pivot, spid);        
}
