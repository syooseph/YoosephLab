#include "rjoinerSA.h"

ReadJoinerSA::ReadJoinerSA()
{
    lend_gsa = NULL;
    rend_gsa = NULL;
    lend_path_strs = NULL;
    rend_path_strs = NULL;
    lend_path_lens = NULL;
    rend_path_lens = NULL;

    if ( Param::pair_flag && Param::rjoin_pe_support )
    bridge_log.setPESupport(true);
    max_read_len = 0;
}

ReadJoinerSA::~ReadJoinerSA()
{
    purgeTemp();
    // if ( lend_gsa != NULL ) delete lend_gsa;
    // if ( rend_gsa != NULL ) delete rend_gsa;
}

void ReadJoinerSA::initMaps()
{
    getMaxReadLength();
    extractEndStrings();
    makeSuffixArrays();
    extractBridgeReads();
}

void ReadJoinerSA::getMaxReadLength()
{
    double t0 = mytime();
    max_read_len = 0;
    for ( int i = 0; i < nreads; i++ ) {
        int l = strlen(seqs[i]);
        if ( l > max_read_len ) 
            max_read_len = l;
    }
    if ( Param::verbose >= 1 ) 
        printf( "Max. read length:%d\ttime:%.4f\n", max_read_len, mytime()-t0 );
}

// pure virutal function from Connectr class
bool ReadJoinerSA::latchOverlapPaths( PathId &spid, int direction )
{
    return true;
}

// pure virutal function from Connectr class
void ReadJoinerSA::updateIndex( PathId spid,
                                 std::string &old_pivot,
                                 std::string &new_pivot,
                                 int direction)
{

}


void ReadJoinerSA::extractEndStrings()
{
    if ( Param::ncpus == 1 ) {
        extractOneEndStrings(LEFT);
        extractOneEndStrings(RIGHT);
    } else {
        std::thread t1(&ReadJoinerSA::extractOneEndStrings, this, LEFT);
        std::thread t2(&ReadJoinerSA::extractOneEndStrings, this, RIGHT);
        t1.join();
        t2.join();
    }
}

void ReadJoinerSA::extractOneEndStrings( int direction ) 
{
    double t0 = mytime();

    size_t n = extenders.size();

    if ( direction == LEFT ) {
        lend_path_strs = new char*[n];
        lend_path_lens = new size_t[n];
    }
    else {
        rend_path_strs = new char*[n];
        rend_path_lens = new size_t[n];
    }

    PathId nid = 0;
    for ( auto item : extenders ) {
        PathId pid = item.first;
        std::string str = item.second.sequence;

        assert( str.size() > 0 );
        //if ( str.size() < max_read_len ) continue;

        size_t len = (int)str.size() >= max_read_len ?
            max_read_len : str.size();
        
        std::string end_str = direction == LEFT ?
            str.substr(0,len) :
            str.substr(str.size()-len,len);

        if ( direction == LEFT ) {
            lend_path_strs[nid] = new char[end_str.size()+1];
            strcpy( lend_path_strs[nid], end_str.c_str() );
            lend_path_lens[nid] = len;
        } else {
            rend_path_strs[nid] = new char[end_str.size()+1];
            strcpy( rend_path_strs[nid], end_str.c_str() );
            rend_path_lens[nid] = len;
        }

        // need only once
        if ( direction == LEFT ) 
            gsa_ids.insert( std::pair<PathId, PathId>( nid, pid ) );
        nid++;
    }
    
    if ( Param::verbose >= 1 ) 
        printf( "Path end strings (direction:%d) extracted:%.4f\n", direction, mytime()-t0 );
}

void ReadJoinerSA::makeSuffixArrays()
{
    double t0 = mytime();

    assert( Param::ncpus >= 1 );

    if ( Param::ncpus == 1 ) {
        makeOneSuffixArray( LEFT );
        makeOneSuffixArray( RIGHT );
    } else {
        std::thread t1(&ReadJoinerSA::makeOneSuffixArray, this, LEFT);
        std::thread t2(&ReadJoinerSA::makeOneSuffixArray, this, RIGHT);
        t1.join();
        t2.join();
    }

    if ( Param::verbose >= 1 ) 
        printf( "Suffix arrays created:%.4f\n", mytime()-t0 );
}

void ReadJoinerSA::makeOneSuffixArray( int direction )
{
    double t0 = mytime();

    if ( direction == LEFT ) 
        lend_gsa = new GSA( lend_path_strs, gsa_ids.size(), true );
    else
        rend_gsa = new GSA( rend_path_strs, gsa_ids.size(), true );
 
    // for ( auto item : sub_strs ) 
    //     delete[] item;
    
    if ( Param::verbose >= 1 ) 
        printf( "Suffix array (direction:%d) created:%.4f\n", direction, mytime()-t0 );
}

int ReadJoinerSA::countUnrecruitedReads()
{
    assert( preads != NULL );
    if ( Param::verbose >= 1 ) 
        std::cout << "# reads:" << nreads << "\n";
    int nquery = 0;
    for ( int i = 0; i < nreads; i++ ) {
        if ( preads[i] != NOT_PATH ) continue;
        if ( preads[i] == BAD_READ ) continue;
        if ( (int)strlen( seqs[i] ) < Param::bridge_overlap ) continue;
        nquery++;
    }
    if ( Param::verbose >= 1 ) 
        std::cout << "# unrecruited:" << nquery << "\n";
    return nquery;
}

void ReadJoinerSA::extractBridgeReads()
{
    std::cout << "Searching all bridging reads ...\n";
    Param::ncpus < 2 ? __extractBridgeReadsSP() : __extractBridgeReadsMP();
}

void ReadJoinerSA::__extractBridgeReadsMP()
{
    double t0 = mytime();

    int nquery = countUnrecruitedReads();
    Progress prog( 1.0, 0, nquery, mytime() );

    int l = Param::bridge_overlap;
    int i;

    std::vector< ReadToPathsMap > LendReadPaths( Param::ncpus, ReadToPathsMap() );
    std::vector< ReadToPathsMap > RendReadPaths( Param::ncpus, ReadToPathsMap() );
    
    std::vector< ReadMatchPossMap > LendPathReads( Param::ncpus, ReadMatchPossMap() );
    std::vector< ReadMatchPossMap > RendPathReads( Param::ncpus, ReadMatchPossMap() );

    t_suffix = t_trim = t_align = t_save = 0.0;
    t_trim_right = t_trim_left = t_trim_erase = 0.0;
    t_align_check = t_align_insert = t_align_rlen =  t_align_getsa = t_align_plen = 0.0;
    n_align_total = n_align_check = n_align_valid = 0;

    double prev = prog.ratio;

#pragma omp parallel for schedule(dynamic, 1) private(i) num_threads(Param::ncpus)
    for ( i = 0; i < nreads; i++ ) {

        int tid = omp_get_thread_num(); 

        if ( preads[i] != NOT_PATH ) continue;
        if ( preads[i] == BAD_READ ) continue;
        if ( (int)strlen( seqs[i] ) < Param::bridge_overlap ) continue;

#pragma omp critical 
{
        prog.count++;
        prog.showProgress();
        if ( prog.ratio > prev && Param::verbose >= 1 ) {
            printf("- search:%.4f, align:%.4f (rlen:%.4f, getsfa:%.4f, get-plen:%.4f, check:%.4f, insert:%.4f), save:%.4f\n",
                   t_suffix, t_align, t_align_rlen, t_align_getsa, t_align_plen, t_align_check, t_align_insert, t_save );
            printf("- #total:%zu, #check:%zu, #valid:%zu\n", n_align_total, n_align_check, n_align_valid );
            prev = prog.ratio;
        }
 }

        int rlen = strlen(seqs[i]);

        double tic = mytime();
        BoundType rsuffix = rend_gsa->search((SfaChar*)&seqs[i][0], l);
        BoundType lsuffix = lend_gsa->search((SfaChar*)&seqs[i][rlen-l], l);
        //#pragma omp critical
        t_suffix += (mytime()-tic);

        if ( rsuffix.second < rsuffix.first ) {
            if ( Param::verbose >= 3 ) std::cout << "rid:" << i << "\tskipped (rsuffix)\n";
            continue;
        }
        if ( lsuffix.second < lsuffix.first ) {
            if ( Param::verbose >= 3 ) std::cout << "rid:" << i << "\tskipped (lsuffix)\n";
            continue;
        }
        if ( Param::verbose >= 3 ) {
            char prefix[l+1];
            strncpy(prefix, &seqs[i][0], l);
            prefix[l] = '\0';
            std::cout << i << "\t" << seqs[i] << "\n";
            std::cout << "prefix:" << prefix << "\tsuffix:" << &seqs[i][rlen-l] << "\n";
            std::cout << "lsuffix:" << lsuffix.first << ", " << lsuffix.second << "\t" << lsuffix.second-lsuffix.first+1 << "\n";
            std::cout << "rsuffix:" << rsuffix.first << ", " << rsuffix.second << "\t" << rsuffix.second-rsuffix.first+1 << "\n";
        }

        // tic = mytime();
        // dropCommonPaths( lsupports, rsupports, lsuffix, rsuffix );
        // //#pragma omp critical
        // t_trim += (mytime()-tic);
        // if ( lsupports.size() == 0 || rsupports.size() == 0 ) {
        //     if ( Param::verbose ) std::cout << "rid:" << i << "\tskipped (low supports)\n";
        //     continue;
        // }

        PathPosInfoMap rgood_paths, lgood_paths;
        tic = mytime();
        extractValidPaths( i, lgood_paths, rgood_paths, lsuffix, rsuffix );
        //#pragma omp critical
        t_align += (mytime()-tic);
        // if ( rgood_paths.size() == 0 || lgood_paths.size() == 0 ) {
        //     if ( Param::verbose ) std::cout << "rid:" << i << "\tskipped (low valids)\n";
        //     continue;
        // }


        tic = mytime();
        for ( auto item : lgood_paths ) {
            PathId  sid = item.first;
            assert( gsa_ids.find(sid) != gsa_ids.end() );
            PathId pid = gsa_ids[sid];
            PosInfo inf = item.second;
            LendReadPaths[tid][(ReadId)i].insert(pid);
            ReadMatchPos rmp( (ReadId)i, inf.rpos, inf.ppos, inf.mlen );
            LendPathReads[tid][pid].push_back( rmp );
            
            if ( Param::verbose >= 3 ) {
                printf("lend-rid:%d, pid:%d, rpos:%d, ppos:%d, mlen:%d\n", i, pid, inf.rpos, inf.ppos, inf.mlen );
                char np[inf.mlen+1];
                strncpy( np, &extenders[pid].sequence[inf.ppos], inf.mlen );
                np[inf.mlen] = '\0';
                
                char nr[inf.mlen+1];
                strncpy( nr, &seqs[i][inf.rpos], inf.mlen );
                nr[inf.mlen] = '\0';
                std::cout << "npath:" << np << "\n";
                std::cout << "nread:" << nr << "\n";
            }
        }

        for ( auto item : rgood_paths ) {
            PathId  sid = item.first;
            assert( gsa_ids.find(sid) != gsa_ids.end() );
            PathId pid = gsa_ids[sid];
            //PathId  pid = item.first;
            PosInfo inf = item.second;
            RendReadPaths[tid][(ReadId)i].insert(pid);
            int npos = extenders[pid].sequence.size() - inf.mlen;
            assert(npos >= 0 );
            ReadMatchPos rmp( (ReadId)i, inf.rpos, npos, inf.mlen );
            RendPathReads[tid][pid].push_back( rmp );

            if ( Param::verbose >= 3 ) {
                printf("rend-rid:%d, pid:%d, rpos:%d, ppos:%d, mlen:%d\n", i, pid, inf.rpos, npos, inf.mlen );
                char np[inf.mlen+1];
                strncpy( np, &extenders[pid].sequence[npos], inf.mlen );
                np[inf.mlen] = '\0';
                
                char nr[inf.mlen+1];
                strncpy( nr, &seqs[i][inf.rpos], inf.mlen );
                nr[inf.mlen] = '\0';
                std::cout << "npath:" << np << '\n';
                std::cout << "nread:" << nr << "\n";
            }
        }
        //#pragma omp critical
        t_save += (mytime()-tic);
    }

    if ( Param::verbose >= 1 ) {
        printf( "Bridging reads extracted:%.4f\n", mytime()-t0 );
    }

    t0 = mytime();
    for ( int i = 0; i < Param::ncpus; i++ ) {
        for ( auto item : LendReadPaths[i] )
            lend_read_paths.insert( std::pair<ReadId, PathIdSet>( item ) );
        LendReadPaths[i].clear();
        
        for ( auto item : RendReadPaths[i] )
            rend_read_paths.insert( std::pair<ReadId, PathIdSet>( item ) );
        RendReadPaths[i].clear();

        for ( auto item : LendPathReads[i] )
            for ( auto jtem : item.second ) 
                lend_path_reads[item.first].push_back( jtem );
        LendPathReads[i].clear();

        for ( auto item : RendPathReads[i] )
            for ( auto jtem : item.second ) 
                rend_path_reads[item.first].push_back( jtem );
        RendPathReads[i].clear();
    }

    if ( Param::verbose >= 1 ) {
        std::cout << "Lend reads:" << lend_read_paths.size() << "\n";
        std::cout << "Rend reads:" << rend_read_paths.size() << "\n";
        std::cout << "Lend paths:" << lend_path_reads.size() << "\n";
        std::cout << "Rend paths:" << rend_path_reads.size() << "\n";
    }

    if ( Param::verbose >= 1 ) 
        printf( "Bridging reads combined:%.4f\n", mytime()-t0 );

    t0 = mytime();
    ReadMatchPossMap::iterator it;
    for ( it = rend_path_reads.begin(); it != rend_path_reads.end(); )
        (int)it->second.size() < Param::min_bridge_reads ? rend_path_reads.erase(it++) : ++it;
    
    for ( it = lend_path_reads.begin(); it != lend_path_reads.end(); )
        (int)it->second.size() < Param::min_bridge_reads ? lend_path_reads.erase(it++) : ++it;

    if ( Param::verbose >= 1 ) {
        printf("Weak support paths dropped:%.4f\n", mytime()-t0);
        std::cout << "Lend paths:" << lend_path_reads.size() << "\n";
        std::cout << "Rend paths:" << rend_path_reads.size() << "\n";
    }
}


void ReadJoinerSA::__extractBridgeReadsSP()
{
    double t0 = mytime();

    int nquery = countUnrecruitedReads();
    Progress prog( 1.0, 0, nquery, mytime() );

    int l = Param::bridge_overlap;
    int i;

    t_suffix = t_trim = t_align = t_save = 0.0;
    t_trim_right = t_trim_left = t_trim_erase = 0.0;
    t_align_check = t_align_insert = t_align_rlen =  t_align_getsa = t_align_plen = 0.0;
    n_align_total = n_align_check = n_align_valid = 0;

    double prev = prog.ratio;

    for ( i = 0; i < nreads; i++ ) {

        if ( preads[i] != NOT_PATH ) continue;
        if ( preads[i] == BAD_READ ) continue;
        if ( (int)strlen( seqs[i] ) < Param::bridge_overlap ) continue;

        prog.count++;
        prog.showProgress();
        if ( prog.ratio > prev && Param::verbose >= 1 ) {
            printf("- search:%.4f, align:%.4f (rlen:%.4f, getsfa:%.4f, get-plen:%.4f, check:%.4f, insert:%.4f), save:%.4f\n",
                   t_suffix, t_align, t_align_rlen, t_align_getsa, t_align_plen, t_align_check, t_align_insert, t_save );
            printf("- #total:%zu, #check:%zu, #valid:%zu\n", n_align_total, n_align_check, n_align_valid );
            prev = prog.ratio;
        }

        int rlen = strlen(seqs[i]);

        double tic = mytime();
        BoundType rsuffix = rend_gsa->search((SfaChar*)&seqs[i][0], l);
        BoundType lsuffix = lend_gsa->search((SfaChar*)&seqs[i][rlen-l], l);
        t_suffix += (mytime()-tic);

        if ( rsuffix.second < rsuffix.first ) {
            if ( Param::verbose >= 3 ) std::cout << "rid:" << i << "\tskipped (rsuffix)\n";
            continue;
        }
        if ( lsuffix.second < lsuffix.first ) {
            if ( Param::verbose >= 3 ) std::cout << "rid:" << i << "\tskipped (lsuffix)\n";
            continue;
        }
        if ( Param::verbose >= 3 ) {
            char prefix[l+1];
            strncpy(prefix, &seqs[i][0], l);
            prefix[l] = '\0';
            std::cout << i << "\t" << seqs[i] << "\n";
            std::cout << "prefix:" << prefix << "\tsuffix:" << &seqs[i][rlen-l] << "\n";
            std::cout << "lsuffix:" << lsuffix.first << ", " << lsuffix.second << "\t" << lsuffix.second-lsuffix.first+1 << "\n";
            std::cout << "rsuffix:" << rsuffix.first << ", " << rsuffix.second << "\t" << rsuffix.second-rsuffix.first+1 << "\n";
        }

        PathPosInfoMap rgood_paths, lgood_paths;
        tic = mytime();
        extractValidPaths( i, lgood_paths, rgood_paths, lsuffix, rsuffix );
        t_align += (mytime()-tic);

        tic = mytime();
        for ( auto item : lgood_paths ) {
            PathId  sid = item.first;
            assert( gsa_ids.find(sid) != gsa_ids.end() );
            PathId pid = gsa_ids[sid];
            PosInfo inf = item.second;
            lend_read_paths[(ReadId)i].insert(pid);
            ReadMatchPos rmp( (ReadId)i, inf.rpos, inf.ppos, inf.mlen );
            lend_path_reads[pid].push_back( rmp );            
            if ( Param::verbose >= 3 ) {
                printf("lend-rid:%d, pid:%d, rpos:%d, ppos:%d, mlen:%d\n", i, pid, inf.rpos, inf.ppos, inf.mlen );
                char np[inf.mlen+1];
                strncpy( np, &extenders[pid].sequence[inf.ppos], inf.mlen );
                np[inf.mlen] = '\0';
                
                char nr[inf.mlen+1];
                strncpy( nr, &seqs[i][inf.rpos], inf.mlen );
                nr[inf.mlen] = '\0';
                std::cout << "npath:" << np << "\n";
                std::cout << "nread:" << nr << "\n";
            }
        }

        for ( auto item : rgood_paths ) {
            PathId  sid = item.first;
            assert( gsa_ids.find(sid) != gsa_ids.end() );
            PathId pid = gsa_ids[sid];
            //PathId  pid = item.first;
            PosInfo inf = item.second;
            rend_read_paths[(ReadId)i].insert(pid);
            int npos = extenders[pid].sequence.size() - inf.mlen;
            assert(npos >= 0 );
            ReadMatchPos rmp( (ReadId)i, inf.rpos, npos, inf.mlen );
            rend_path_reads[pid].push_back( rmp );
            if ( Param::verbose >= 3 ) {
                printf("rend-rid:%d, pid:%d, rpos:%d, ppos:%d, mlen:%d\n", i, pid, inf.rpos, npos, inf.mlen );
                char np[inf.mlen+1];
                strncpy( np, &extenders[pid].sequence[npos], inf.mlen );
                np[inf.mlen] = '\0';
                
                char nr[inf.mlen+1];
                strncpy( nr, &seqs[i][inf.rpos], inf.mlen );
                nr[inf.mlen] = '\0';
                std::cout << "npath:" << np << '\n';
                std::cout << "nread:" << nr << "\n";
            }
        }

        t_save += (mytime()-tic);
    }

    if ( Param::verbose >= 1 ) {
        printf( "Bridging reads extracted:%.4f\n", mytime()-t0 );
    }

    if ( Param::verbose >= 1 ) {
        std::cout << "Lend reads:" << lend_read_paths.size() << "\n";
        std::cout << "Rend reads:" << rend_read_paths.size() << "\n";
        std::cout << "Lend paths:" << lend_path_reads.size() << "\n";
        std::cout << "Rend paths:" << rend_path_reads.size() << "\n";
    }

    t0 = mytime();
    ReadMatchPossMap::iterator it;
    for ( it = rend_path_reads.begin(); it != rend_path_reads.end(); )
        (int)it->second.size() < Param::min_bridge_reads ? rend_path_reads.erase(it++) : ++it;
    
    for ( it = lend_path_reads.begin(); it != lend_path_reads.end(); )
        (int)it->second.size() < Param::min_bridge_reads ? lend_path_reads.erase(it++) : ++it;

    if ( Param::verbose >=1 ) {
        printf("Weak support paths dropped:%.4f\n", mytime()-t0);
        std::cout << "Lend paths:" << lend_path_reads.size() << "\n";
        std::cout << "Rend paths:" << rend_path_reads.size() << "\n";
    }
}

void ReadJoinerSA::dropCommonPaths( PathPossMap &lsupports, 
                                    PathPossMap &rsupports, 
                                    const BoundType &lsuffix,
                                    const BoundType &rsuffix )
{
    double t0 = mytime();
    for ( int j = rsuffix.first; j <= rsuffix.second; j++ ) {
        GsaType sa = rend_gsa->getAt(j);
        PathId  sid = sa.doc;
        assert( gsa_ids.find(sid) != gsa_ids.end());
        int     pos = (int)sa.pos;

        if ( Param::verbose ) {
            std::cout << "rend-sid:" << sid << "\tpid:" << gsa_ids[sid] << "\tpos:" << pos << "\trend-seq:" << &rend_path_strs[sid][pos] << "\n";
            std::cout << "orig path:" << extenders[gsa_ids[sid]].sequence << "\n";
        }

        double tic = mytime();
        rsupports[sid].push_back(pos);
        //#pragma omp critical
        t_trim_insert += (mytime()-tic);
    }
    //#pragma omp critical
    t_trim_right += (mytime()-t0);
    
    if ( Param::verbose ) std::cout << "# right:" << rsupports.size() << "\n";

    t0 = mytime();
    for ( int j = lsuffix.first; j <= lsuffix.second; j++ ) {
        GsaType sa = lend_gsa->getAt(j);
        PathId  sid = sa.doc;
        assert( gsa_ids.find(sid) != gsa_ids.end());
        int     pos = (int)sa.pos;

        //if ( rsupports.find(sid) != rsupports.end() ) continue;
        if ( Param::verbose ) {
            std::cout << "lend-sid:" << sid << "\tpid:" << gsa_ids[sid] << "\tpos:" << pos << "\tlend-seq:" << &lend_path_strs[sid][pos] << "\n";
            std::cout << "orig path:" << extenders[gsa_ids[sid]].sequence << "\n";
        }

        double tic = mytime();
        lsupports[sid].push_back(pos);
        //#pragma omp critical
        t_trim_insert += (mytime()-tic);

    }
    //#pragma omp critical
    t_trim_left += (mytime()-t0);

    if ( Param::verbose ) std::cout << "# left:" << lsupports.size() << "\n";

//     t0 = mytime();
//     std::tr1::unordered_map< PathId, IntArray>::iterator rt;
//     for ( rt = rsupports.begin(); rt != rsupports.end(); ) {
//         if ( lsupports.find(rt->first) != lsupports.end() )
//             rt = rsupports.erase(rt++);
//         else ++rt;
//     }
// #pragma omp critical
//     t_trim_erase += (mytime()-t0);


    if ( Param::verbose ) {
        std::cout << "lsupports:" << lsupports.size() << "\n";
        std::cout << "rsupports:" << rsupports.size() << "\n";
    }
}

void ReadJoinerSA::extractValidPaths( int i,
                                      PathPosInfoMap &lgood_paths,
                                      PathPosInfoMap &rgood_paths, 
                                      const BoundType &lsuffix,
                                      const BoundType &rsuffix ) 
{
    double t0 = mytime();
    int rlen = strlen(seqs[i]);
    t_align_rlen += (mytime()-t0);

    int l = Param::bridge_overlap;
    n_align_total += (lsuffix.second-lsuffix.first+1);
    n_align_total += (rsuffix.second-rsuffix.first+1);

    // check string match
    // .....ooooo??        path
    //      ooooo??....... read
    for ( int j = rsuffix.first; j <= rsuffix.second; j++ ) {
        double tic = mytime();
        GsaType sa = rend_gsa->getAt(j);
        t_align_getsa += (mytime()-tic);

        PathId  sid = sa.doc;
        if ( rgood_paths.find(sid) != rgood_paths.end() ) continue;
        //assert( gsa_ids.find(sid) != gsa_ids.end());
        //PathId pid = gsa_ids[sid];
        int    pos = (int)sa.pos;
        
        tic = mytime();
        //int plen = strlen(rend_path_strs[sid]);
        int plen = rend_path_lens[sid];
        t_align_plen += (mytime()-tic);

        int roff = rlen - l;
        int poff = plen - (pos+l);
        if ( Param::verbose >= 3 ) 
            printf("rend-sid:%d, pid:%d, plen:%d, roff:%d, poff:%d\n", sid, gsa_ids[sid], plen, roff, poff );
        if ( roff < poff ) continue;
        
        int ps = pos+l;
        int rs = l;
        
        if ( Param::verbose >= 3 ) {
            printf("rs:%d, ps:%d, rl:%d, pl:%d, len:%d\n", rs, ps, rlen, plen, poff );
            char rsub[poff+1];
            strncpy( rsub, &seqs[i][rs], poff );
            rsub[poff] = '\0';
            char psub[poff+1];
            strncpy( psub, &rend_path_strs[sid][ps], poff );
            psub[poff] = '\0';
            std::cout << "rsub:" << rsub << "\tpsub:" << psub << "\n";
        }

        tic = mytime();
        bool good = sameBases( seqs[i], rend_path_strs[sid], rs, ps, rlen, plen, poff );
        //#pragma omp critical
        t_align_check += (mytime()-tic);
        n_align_check++;
        if ( Param::verbose >= 3 ) std::cout << "same base (rend)?:" << good << "\n";
        if ( !good ) continue;
        
        n_align_valid++;
        tic = mytime();                  
        PosInfo p(0, pos, l+poff);
        if ( Param::verbose >= 3 ) 
            printf("rs:%d, ps:%d, l:%d\n", 0, pos, l+poff);
        
        rgood_paths.insert( std::pair<PathId, PosInfo>( sid, p ) );
        //#pragma omp critical
            t_align_insert += (mytime()-tic);
    }
    
    // check string match
    //      ??ooooo....... path
    // .....??ooooo        read
    for ( int j = lsuffix.first; j <= lsuffix.second; j++ ) {
        double tic = mytime();
        GsaType sa = lend_gsa->getAt(j);
        t_align_getsa += (mytime()-tic);

        PathId  sid = sa.doc;
        if ( lgood_paths.find(sid) != lgood_paths.end() ) continue;
        //assert( gsa_ids.find(sid) != gsa_ids.end());
        //PathId pid = gsa_ids[sid];
        int     pos = (int)sa.pos;
        

        tic = mytime();
        //int plen = strlen(lend_path_strs[sid]);
        int plen = lend_path_lens[sid];
        t_align_plen += (mytime()-tic);
        // int rlen = strlen(seqs[i]);
        int roff = rlen - l;
        int poff = pos;
        if ( Param::verbose >= 3 ) 
            printf("lend-sid:%d, pid:%d, plen:%d, roff:%d, poff:%d\n", sid, gsa_ids[sid], plen, roff, poff );
        if ( roff < poff ) continue;
        
        int rs = rlen-l-poff;
        int ps = 0;

        if ( Param::verbose >= 3 ) {
            printf("rs:%d, ps:%d, rl:%d, pl:%d, len:%d\n", rs, ps, rlen, plen, poff );
            char rsub[poff+1];
            strncpy( rsub, &seqs[i][rs], poff );
            rsub[poff] = '\0';
            char psub[poff+1];
            strncpy( psub, &lend_path_strs[sid][ps], poff );
            psub[poff] = '\0';
            std::cout << "rsub:" << rsub << "\tpsub:" << psub << "\n";
        }
        
        int ors = rs, ops = ps;
        tic = mytime();
        bool good = sameBases( seqs[i], lend_path_strs[sid], rs, ps, rlen, plen, poff );
        //#pragma omp critical
        t_align_check += (mytime()-tic);
        n_align_check++;
        if ( Param::verbose >= 3 ) std::cout << "same base (lend)?:" << good << "\n";
        if ( !good ) continue;
        n_align_valid++;

        tic = mytime();
        PosInfo p(ors, ops, l+poff);
        if ( Param::verbose >= 3 ) 
            printf("rs:%d, ps:%d, l:%d\n", ors, ops, l+poff);
        
        lgood_paths.insert( std::pair<PathId, PosInfo>( sid, p ) );
        //#pragma omp critical
            t_align_insert += (mytime()-tic);;
    }

    
    if ( Param::verbose >= 3 ) {
        std::cout << "lvalid:" << lgood_paths.size() << "\n";
        std::cout << "rvalid:" << rgood_paths.size() << "\n";
    }

}

bool ReadJoinerSA::sameBases( const char *rstr,
                              const char *pstr,
                              int &rs,
                              int &ps,
                              const int &rlen,
                              const int &plen,
                              int &len )
{
    // int l1 = strlen(str1);
    // int l2 = strlen(str2);
    for ( int i = 0; i<len; i++,rs++,ps++) {
        //if ( s1 >= l1 ) return false;
        //if ( s2 >= l2 ) return false;
        assert(rs<rlen);
        assert(ps<plen);
        if ( rstr[rs] != pstr[ps] ) 
            return false;
    }

    return true;
}


//void ReadJoinerSA::stitchPaths()
void ReadJoinerSA::joinPaths()
{
    std::cout << "\nExtending paths ...\n";
    size_t npath  = extenders.size();
    progress = Progress( 1.0, 0, npath, mytime() );

    double t0 = mytime();

    PathLengthMap plen_map = getPathLengths();
    for ( auto it = plen_map.rbegin(); it != plen_map.rend(); ++it ) {
        /* Print progress status at every percent of path proceeded */
        double prev = progress.ratio;
        progress.count++;
        progress.showProgress();
        bridge_log.et_total = (mytime()-t0);
        if ( Param::verbose >= 1 && progress.ratio > prev ) 
            if ( progress.ratio > prev ) bridge_log.printSummary();

        PathId sbjct_pid = it->second;
        if ( merged_paths.find( sbjct_pid ) != merged_paths.end() ) continue;

        stitchPath( sbjct_pid, LEFT );
        stitchPath( sbjct_pid, RIGHT );
    }
    
    bridge_log.et_total = mytime()-t0;
    
    std::cout << "# Merged paths:" << merged_paths.size() << "\n";
    std::cout << "# Bridging reads:" << getAddedReadsCount() << "\n";
    
    //if ( Param::verbose ) printAddedReads();
    if ( Param::verbose >= 3 || Param::debug_flag ) printAddedReads();
    
}

void ReadJoinerSA::stitchPath( PathId spid, int direction )
{
    if ( latchReadBridgingPaths( spid, direction ) )
        stitchPath( spid, direction );
    else return;
}


bool ReadJoinerSA::latchReadBridgingPaths( PathId spid, int direction )
{
    double t0 = mytime();

    std::multimap<size_t, PathId> match_counts;
    findCandidates( match_counts, spid, direction );
    bridge_log.et_candidate += (mytime()-t0);
    bridge_log.ct_candidates += match_counts.size();
    if ( match_counts.size() == 0 ) return false;

    if ( match_counts.size() > bridge_log.ct_cand_max ) bridge_log.ct_cand_max = match_counts.size();
    if ( match_counts.size() < bridge_log.ct_cand_min ) bridge_log.ct_cand_min = match_counts.size();
    bridge_log.ct_iter++;

    double tic = mytime();
    /* Find reads in library size length */
    std::tr1::unordered_map<ReadId, bool> pivot_reads;
    if ( Param::pair_flag && Param::rjoin_pe_support ) {
        std::string pivot = extenders[spid].sequence;
        std::string sub_pivot = getStringInRange(pivot, direction, true );            
        findReads( sub_pivot, pivot_reads );
    }
    bridge_log.et_read_pair += (mytime()-tic);
    if ( Param::verbose ) 
        printf("# pivot reads in library range:%zu (%.4f sec)\n", pivot_reads.size(), mytime()-tic );

    tic = mytime();
    bool success = false;
    size_t iter = 0;
    size_t max_iter = 10;
    std::multimap<size_t, PathId>::reverse_iterator rt;
    for ( rt = match_counts.rbegin(); rt != match_counts.rend(); ++rt ) {
        ++iter;
        if ( iter > max_iter ) break;
        if ( (int)rt->first < Param::min_bridge_reads ) break;
     
        
        ReadMatchPossMap  *pivot_path_reads = ( direction == LEFT )  ?
            &lend_path_reads : &rend_path_reads; 
        
        ReadMatchPossMap::iterator pt = pivot_path_reads->find(spid);
        assert( pt != pivot_path_reads->end() );
        ReadMatchPosList *pivot_poss = &(*pivot_path_reads)[spid];

        PathId mpid = rt->second;
        bridge_log.ct_trial++;
        success = tryLatch( spid, mpid, pivot_poss, pivot_reads, direction );
        success ? bridge_log.ct_success++ : bridge_log.ct_fail++;
        if ( iter == 1 ) {
            success ? bridge_log.ct_succ_first++ : bridge_log.ct_fail_first++;
        }
        if ( success ) break;
    }

    if ( Param::verbose ) 
        printf("iter:%zu\tsuccess:%d\ttime:%.4f\n", iter, success, mytime()-tic);

    return success;
}

void ReadJoinerSA::findCandidates( std::multimap<size_t, PathId> &count_map,
                                PathId spid,
                                int direction )
{
    ReadMatchPossMap *pivot_reads = ( direction == LEFT ) ? &lend_path_reads : &rend_path_reads;
    ReadMatchPossMap *match_reads = ( direction == LEFT ) ? &rend_path_reads : &lend_path_reads;

    if ( pivot_reads->find( spid ) == pivot_reads->end() ) return;
    if ( (int)(*pivot_reads)[spid].size() < Param::min_bridge_reads ) return;

    ReadToPathsMap *match_paths = ( direction == LEFT ) ? &rend_read_paths : &lend_read_paths;


    std::tr1::unordered_map<PathId, size_t> match_counts;
    //GsaTypeList valids;
    assert( pivot_reads->find(spid) != pivot_reads->end() );
    for ( auto it = (*pivot_reads)[spid].begin(); it != (*pivot_reads)[spid].end(); ++it ) {
        ReadId rid = it->read;
        ReadToPathsMap::iterator jt = match_paths->find( rid );
        if ( jt == match_paths->end() ) continue;
        for ( auto kt = jt->second.begin(); kt != jt->second.end(); ++kt ) {
            if ( *kt == spid ) continue;
            if ( merged_paths.find(*kt) != merged_paths.end() ) continue;
            // already consumed
            if ( match_reads->find(*kt) == match_reads->end() ) continue;
            // should be ok 
            //match_counts[*kt]++;
            if ( match_counts.find( *kt ) == match_counts.end() ) 
                match_counts.insert( std::pair<PathId, size_t>( *kt, 0 ) );
            match_counts[*kt]++;
        }
    }

    if ( match_counts.size() == 0 ) return;

    if ( Param::verbose ) {
        std::cout << "\nPath:" << spid 
                  << "\t#Candidates:" << match_counts.size() 
                  << "\tdirection:" << direction << "\n";
        for ( auto item : match_counts ) {
            std::cout << "spid:" << spid << "\tmpid:" << item.first << "\tsupport-read-count:" << item.second << "\n";
        }
    }

    double tic = mytime();
    count_map = util::sortByValue<PathId, size_t>( match_counts );
    bridge_log.et_sort += (mytime()-tic);
}


bool ReadJoinerSA::tryLatch( PathId spid,
                             PathId mpid,
                             const ReadMatchPosList *pivot_poss,
                             std::tr1::unordered_map<ReadId, bool> &pivot_reads,
                             int direction )
{
    bool good = false;

    std::string pivot_str = extenders[spid].sequence;
    std::string match_str = extenders[mpid].sequence;

    if ( Param::verbose ) {
        std::cout << "Checking latchability - direction:" << direction << "\n";
        std::cout << "Pivot:" << spid << "\t" << pivot_str << "\n";
        std::cout << "Match:" << mpid << "\t" << match_str << "\n";
    }

    ReadMatchPossMap *match_path_reads = ( direction == LEFT ) ?
        &rend_path_reads : &lend_path_reads;
    
    ReadMatchPossMap::iterator it = match_path_reads->find(mpid);
    if ( it == match_path_reads->end()) {
        std::cerr << "mpid:" << mpid << " not found in match_path_reads\n";
    }
    assert( it != match_path_reads->end());
    ReadMatchPosList *match_poss = &(*match_path_reads)[mpid];

    bridge_log.ct_read_all += match_poss->size();

    ReadPosInfoMap pivot_rmap, match_rmap;
    ReadFlagMap comm_reads;
    
    double tic = mytime();
    good = extractSharedReads( pivot_poss, match_poss, pivot_rmap, match_rmap, comm_reads, spid, mpid, direction );
    bridge_log.et_reads += (mytime()-tic);
    bridge_log.et_latchable += (mytime()-tic);

    if ( !good ) {
        bridge_log.ct_read_fail++; 
        return false;
    } 

    bridge_log.ct_read_share += comm_reads.size();


    std::string mid_str;
    ReadIdSet   bridges;
    int overlap;
    tic = mytime();
    good = determineBridgingEvidence( bridges, overlap, mid_str, pivot_str, match_str, pivot_rmap, match_rmap, direction);
    bridge_log.et_evidence += (mytime()-tic );
    bridge_log.et_latchable += (mytime()-tic);

    if ( !good ) {
        bridge_log.ct_evidence_fail++;
        return false;
    }

    bridge_log.ct_read_good += bridges.size();

    tic = mytime();


    int nlen =  pivot_str.size()+match_str.size();
    if ( Param::verbose ) std::cout << "New sequence size:" << nlen << "\n";
    if ( overlap > 0 ) nlen -= mid_str.size();
    else nlen += mid_str.size();
    if ( Param::pair_flag && Param::rjoin_pe_support ) {
        if ( nlen >= min_libsize ) {
            double lt0 = mytime();
            size_t npairs = findPairendReads( pivot_reads, pivot_str, match_str, direction );
            bridge_log.et_read_pair += (mytime()-lt0);
            if ( Param::verbose ) 
                printf("# match pairs in library range:%zu (%.4f sec)\n", npairs, mytime()-lt0 );
            
            if ( (int)npairs < Param::min_pair_reads ) {
                if ( Param::verbose ) std::cout << "No read support\n";
                log.ct_fail_pread++;
                return false; //continue;
            }
        } else {
            if ( Param::verbose )
            std::cout << "Too short length for pair-end support check\n";
            log.ct_fail_short++;
            return false;
        }
    }
    
    if ( Param::debug_flag ) {
        printf("\nLatching path\tspid:%d, mpid:%d, direction:%d\n", spid, mpid, direction);
        std::cout << "pivot:" << extenders[spid].sequence << "\n";    
        std::cout << "match:" << extenders[mpid].sequence << "\n";            
    }

    AlignSummary summary;
    adjustAlignPositions( summary, overlap, mid_str, spid, mpid, direction );
    if ( Param::debug_flag )
        printf("lgaps:%d,%d\n", summary.lgap.first, summary.lgap.second);
    
    // NOW LATCH
    updateSequence( spid, mpid, overlap, mid_str, direction );
    if ( Param::debug_flag )
        std::cout << "new pivot:" << extenders[spid].sequence << "\n";    
    
    // call base class function
    updateMembers( spid, mpid, summary, direction );
    updateStop(spid, mpid, direction);
    updateTrim(spid, mpid, direction);

    bridge_log.et_latch += ( mytime()-tic );

    tic = mytime();
    addReads( spid, mpid, bridges, pivot_rmap, match_rmap, summary, direction );
    updateMap( spid, mpid, summary, direction );
    merged_paths.insert(mpid);
    bridge_log.et_update += (mytime()-tic);

    return true;
}

void ReadJoinerSA::addReads( PathId spid,
                             PathId mpid,
                             ReadIdSet &bridges,
                             ReadPosInfoMap & pivot_poss,
                             ReadPosInfoMap & match_poss,
                             AlignSummary &summary,
                             int direction )
{
    
    if ( Param::debug_flag ) {
        size_t n = ( added_reads.find(spid) == added_reads.end() ) ? 0 : added_reads[spid].size();
        std::cout << "Existing reads for spid:" << n << "\n";
    }
    
    //-------------------------------------------------------------
    // In case left extension, update existig added reads for pivot
    //-------------------------------------------------------------
    if ( added_reads.find(spid) != added_reads.end() && direction == LEFT ) {
        if ( Param::verbose ) std::cout << "Existing reads for spid:" << mpid << "\tsize:" << added_reads[spid].size() << "\n";
        ReadEntryList::iterator jt;
        for ( jt = added_reads[spid].begin(); jt != added_reads[spid].end(); ++jt ) {
            int op = jt->ppos;
            jt->ppos += summary.lgap.first;
            if ( Param::debug_flag ) 
                printf("path-pos shifted for rid:%d, old:%d, new:%d\n", jt->read, op, jt->ppos);
        }
    }

    if ( Param::debug_flag ) {
        size_t n = ( added_reads.find(mpid) == added_reads.end() ) ? 0 : added_reads[mpid].size();
        std::cout << "Existing reads for mpid:" << n << "\n";
    }
                                     
    //----------------------------------------------------------------------
    // Add added reads from match. If right extension, update align position
    //----------------------------------------------------------------------
    if ( added_reads.find(mpid) != added_reads.end() ) {
        if ( Param::verbose ) std::cout << "Existing reads for mpid:" << mpid << "\tsize:" << added_reads[mpid].size() << "\n";
        ReadEntryList::iterator jt;
        for ( jt = added_reads[mpid].begin(); jt != added_reads[mpid].end(); ++jt ) {
            if ( direction == RIGHT ) {
                int op = jt->ppos;
                jt->ppos += summary.lgap.second;
                if ( Param::debug_flag ) 
                    printf("path-pos shifted for rid:%d, old:%d, new:%d\n", jt->read, op, jt->ppos);
            }
            added_reads[spid].push_back(*jt);
        }
        added_reads.erase(mpid);
    }


    if ( Param::debug_flag ) 
        std::cout << "New reads for spid:" << bridges.size() << "\n";

    if ( Param::verbose ) std::cout << "New reads for spid:" << spid << "\tsize:" << bridges.size() << "\n";

    //---------------------------
    // Now add new bridging reads
    //---------------------------
    for ( auto rid : bridges ) {
        unsigned rpos, ppos;
        ReadPosInfoMap::iterator pt;
        if ( direction == LEFT ) {
            pt = match_poss.find(rid);
            assert( pt != match_poss.end() );
            rpos = pt->second.rpos;
            ppos = pt->second.ppos;
        } else {
            pt = pivot_poss.find(rid);
            assert( pt != pivot_poss.end() );
            rpos = pt->second.rpos;
            ppos = pt->second.ppos;
        }

        if ( Param::debug_flag || Param::verbose >= 3 ) 
            printf("new read:%d, rpos:%d, ppos:%d\n", rid, rpos, ppos);

        if ( Param::verbose >= 3 ) std::cout << "rid:" << rid << "\trpos:" << rpos << "\tppos:" << ppos << "\n";
        added_reads[spid].push_back( ReadEntry(rid,rpos,ppos) );

        recruited.insert( std::pair<ReadId, PathId>(rid,spid) );

        lend_read_paths.erase(rid);
        rend_read_paths.erase(rid);
    }
}


void ReadJoinerSA::updateMap( PathId spid,
                              PathId mpid,
                              AlignSummary &summary,
                              int direction )
{
    ReadMatchPossMap::iterator it;
    if ( direction == LEFT ) {
        //----------------------------------------------------------------------------
        // Remove left end supporting reads of pivot, which was other latch candiates.
        //----------------------------------------------------------------------------
        it = lend_path_reads.find(spid);
        if ( it != lend_path_reads.end() ) lend_path_reads.erase(it);

        //-------------------------------------------------------------------------------
        // Get left end supporting reads of match.
        // If exist, save them as left end supporting reads of pivot for recursive latch.
        //-------------------------------------------------------------------------------
        it = lend_path_reads.find(mpid);
        if ( it != lend_path_reads.end() ) {
            lend_path_reads.insert( std::pair<PathId,ReadMatchPosList>(spid, it->second) );
            lend_path_reads.erase(it);
        }

        //--------------------------------------------------------------------
        // Pivot path joined a path with leading gaps.
        // So, update aligned position of each read in right supporting reads.
        //--------------------------------------------------------------------
        it = rend_path_reads.find(spid);
        if ( it == rend_path_reads.end() ) return;
        for ( auto &item : it->second ) 
            item.ppos += summary.lgap.first;
    } 
    else {
        //-----------------------------------------------------------------------------
        // Remove right end supporting reads of pivot, which was other latch candiates.
        //-----------------------------------------------------------------------------
        it = rend_path_reads.find(spid);
        if ( it != rend_path_reads.end() ) rend_path_reads.erase(it);

        //--------------------------------------------------------------------------------
        // Get right end supporting reads of match.
        // If exist, save them as right end supporting reads of pivot for recursive latch.
        // Match path was joined into pivot path with leading gaps.
        // So, update aligned position of each read from match path.
        //--------------------------------------------------------------------------------
        it = rend_path_reads.find(mpid);
        if ( it == rend_path_reads.end() ) return;
        for ( auto &item : it->second ) {
            item.ppos += summary.lgap.second;
            rend_path_reads[spid].push_back(item);
        }
        rend_path_reads.erase(it);
    }
}

//====================================================================
// Now, find an evidence to stitch two paths.
// Find a middle string which is either a string to connect two paths 
// or tiny tiny overlap of two paths
//====================================================================
bool ReadJoinerSA::determineBridgingEvidence( ReadIdSet &bridges,
                                            int &overlap,
                                            std::string &mid_str,
                                            std::string &pivot_str,
                                            std::string &match_str,
                                            ReadPosInfoMap &pivot_poss,
                                            ReadPosInfoMap &match_poss,
                                            int direction )
{
    std::unordered_map< std::string, ReadIdSet > ladder_supports;
    std::unordered_map< std::string, int > str_overlaps;

    if ( Param::verbose ) {
        std::cout << "pivot:" << pivot_str << "\n"
        << "match:" << match_str << "\n"
        << "direction:" << direction << "\n";
    }

    // get consistent middle string
    for ( auto item : pivot_poss ) {
        ReadId rid = item.first;
        int pivot_rpos = item.second.rpos;
        int pivot_ppos = item.second.ppos;
        
        ReadPosInfoMap::iterator pt = match_poss.find(rid);
        assert( pt != match_poss.end() );
        int match_rpos = pt->second.rpos;
        int match_ppos = pt->second.ppos;
        
        std::string read_str = seqs[rid];

        
        if ( Param::verbose >= 3 ) {
            std::cout << "rid:" << rid << "\t" << read_str << "\n";
            printf("pivot (rpos:%d, ppos:%d)\n", pivot_rpos, pivot_ppos);
            printf("match (rpos:%d, ppos:%d)\n", match_rpos, match_ppos);
        }

        int o = -1; // no overlap
        //-------------------------------------------
        //                     ooooo..........  pivot
        //           ooooo.....ooooo            read
        // ..........ooooo                      match
        //           | l |     r
        //-------------------------------------------
        if ( direction == LEFT ) {
            assert( match_ppos < (int)match_str.size() );
            int l = (int)match_str.size()-match_ppos;
            int r = pivot_rpos;

            if ( Param::verbose >= 3 ) 
                printf("l:%d, r:%d, direction:%d\n", l,r,direction);

            // no overlap
            if ( r >= l ) {
                mid_str = read_str.substr(l, r-l );
            } else {
                o = abs(l-r);
                if ( Param::verbose >= 3 ) {
                    std::cout << "overlap:" << o << "\n";
                    if ( o >= Param::bridge_overlap ) 
                        std::cerr << "long overlap:" << o << "\n";
                }
                //assert( o < Param::bridge_overlap );
                assert( o <= max_read_len );
                mid_str = pivot_str.substr(0,o);
            }
            if ( Param::verbose >= 3 ) 
                std::cout << "mid-str:" << mid_str << "\n";
        } 
        
        //-------------------------------------------
        // ..........ooooo                      pivot
        //           ooooo.....ooooo            read
        //                     ooooo..........  match
        //           | l |     r
        //-------------------------------------------
        else {
            assert( pivot_ppos < (int)pivot_str.size() );
            int l = (int)pivot_str.size()-pivot_ppos;
            int r = match_rpos;

            if ( Param::verbose >= 3 ) 
                printf("l:%d, r:%d, direction:%d\n", l,r,direction);

            // no overlap
            if ( r >= l ) {
                mid_str = read_str.substr(l, r-l );
            } else {
                o = abs(l-r);
                if ( Param::verbose >= 3 ) {std::cout << "overlap:" << o << "\n";
                    if ( o >= Param::bridge_overlap ) 
                        std::cerr << "long overlap:" << o << "\n";
                }
                //assert( o < Param::bridge_overlap );
                assert( o <= max_read_len );
                mid_str = match_str.substr(0,o);

            }

            if ( Param::verbose >= 3 ) 
                std::cout << "mid-str:" << mid_str << "\n";
        } 

        if ( ladder_supports.find( mid_str ) == ladder_supports.end() )
            ladder_supports.insert( std::pair<std::string, ReadIdSet>( mid_str, ReadIdSet() ) );
        ladder_supports[mid_str].insert(rid);

        str_overlaps[mid_str] = o;
    }
    
    ReadIdSet   max_set;
    std::string max_mid;
    size_t max_num = 0;
    for ( auto item : ladder_supports ) {
        if ( item.second.size() > max_num ) {
            max_num = item.second.size();
            max_set = item.second;
            max_mid = item.first;
        }
    }
    
    if ( (int)max_num < Param::min_bridge_reads ) return false;

    mid_str = max_mid;
    bridges = max_set;
    overlap = str_overlaps[max_mid];

    if ( Param::verbose ) {
        std::cout << "max. mid-str:" << mid_str << "\n";
        std::cout << "#. reads:" << bridges.size() << "\n";
        std::cout << "overlap:" << overlap << "\n";
    }
    return true;    
}


//====================================================================
// Extract bridging reads between two paths
//====================================================================
bool ReadJoinerSA::extractSharedReads( const ReadMatchPosList *pivot_poss,
                                       const ReadMatchPosList *match_poss, 
                                       ReadPosInfoMap &pivot_rmap,
                                       ReadPosInfoMap &match_rmap,
                                       ReadFlagMap &comm_reads,
                                       PathId spid,
                                       PathId mpid,
                                       int direction )
{
    ReadMatchPosList::const_iterator gt;

    //------------------------------------------
    // First, save reads covering pivot sequence
    //------------------------------------------
    double tic = mytime();
    for ( gt = pivot_poss->begin(); gt != pivot_poss->end(); ++gt ) {
        if ( recruited.find(gt->read) != recruited.end() ) continue;
        pivot_rmap.insert( std::pair<ReadId,PosInfo>( gt->read, PosInfo(gt->rpos,gt->ppos,gt->mlen) ) );
    }
    bridge_log.et_read_pivot += (mytime()-tic);
    
    if ( Param::verbose ) std::cout << "pivot-end support reads:" << pivot_rmap.size() << "\n";
    if ( (int)pivot_rmap.size() < Param::min_bridge_reads ) return false;
    
    //------------------------------------------------
    // Find reads covering match sequence.
    // At the same time, find common reads from pivot.
    //------------------------------------------------
    tic = mytime();
    for ( gt = match_poss->begin(); gt != match_poss->end(); ++gt ) {
        if ( recruited.find(gt->read) != recruited.end() ) continue;
        if ( pivot_rmap.find(gt->read) == pivot_rmap.end() ) continue;
        comm_reads.insert( std::pair<ReadId,bool>(gt->read, true) );
        match_rmap.insert( std::pair<ReadId,PosInfo>( gt->read, PosInfo(gt->rpos,gt->ppos,gt->mlen) ) );
    }
    
    if ( Param::verbose ) std::cout << "match-end support reads:" << match_rmap.size() << "\n";
    bridge_log.et_read_match += (mytime()-tic);

    //----------------------------------------
    // We need minimum number of common reads.
    //----------------------------------------
    if ( Param::verbose ) std::cout << "common support reads:" << comm_reads.size() << "\n";
    if ( (int)comm_reads.size() < Param::min_bridge_reads ) return false;

    
    //--------------------------------------------------------------
    // Extract reads that have valid aligned ranges againt two paths
    //--------------------------------------------------------------
    tic = mytime();
    //ReadIdSet bad_reads;
    ReadFlagMap bad_reads;
    ReadPosInfoMap::iterator it, jt;

    for ( it = pivot_rmap.begin(); it != pivot_rmap.end(); ++it) {
        ReadId rid = it->first;
        if ( comm_reads.find(it->first) == comm_reads.end() ) {
            bad_reads.insert( std::pair<PathId, bool>( rid, true) );
            continue;
        }
        int pivot_rpos = it->second.rpos;
        
        jt = match_rmap.find(it->first);
        assert( jt != match_rmap.end() );
        int match_rpos = jt->second.rpos;

        //---------------------------------------------------------
        // Valid
        //                   ooo..........  pivot
        //           ooo.....ooo            read
        // ..........ooo                    match
        //           l       r
        // Invalid
        // We handled this already in long/short overlap extension
        //           ooo.................  pivot
        //           oooo.......           read
        // ...........ooo                  match
        //           rl
        //--------------------------------------------------------
        if ( direction == LEFT ) {
            if ( match_rpos >= pivot_rpos ) {
                bad_reads.insert( std::pair<PathId, bool>( rid, true) );
                continue;
            } 
        } 
        //---------------------------------------------------------
        // Valid
        // ..........ooo                    pivot
        //           ooo.....ooo            read
        //                   ooo..........  match
        //           l       r
        // Invalid
        // We handled this already in long/short overlap extension
        // ...........ooo                  pivot
        //           oooo.......           read
        //           ooo.................  match
        //           rl
        //--------------------------------------------------------
        else {
            if ( match_rpos <= pivot_rpos ) {
                bad_reads.insert( std::pair<PathId, bool>( rid, true) );
                continue;
            }
        }
    }
    bridge_log.et_read_bad += (mytime()-tic);    

    tic = mytime();
    for ( auto item : bad_reads ) {
        pivot_rmap.erase(item.first);
        match_rmap.erase(item.first);
        comm_reads.erase(item.first);
    }
    bridge_log.et_read_trim += (mytime()-tic);

    if ( (int)comm_reads.size() < Param::min_bridge_reads ) return false;

    return true;
}

void ReadJoinerSA::purgeTemp()
{
    if ( lend_gsa != NULL ) delete lend_gsa;
    if ( rend_gsa != NULL ) delete rend_gsa;

    if ( lend_path_lens != NULL ) delete[] lend_path_lens;
    if ( rend_path_lens != NULL ) delete[] rend_path_lens;
    
    if ( lend_path_strs != NULL ) {
        for ( size_t i = 0; i < gsa_ids.size(); i++ ) 
            delete[] lend_path_strs[i];
    }

    if ( rend_path_strs != NULL ) {
        for ( size_t i = 0; i < gsa_ids.size(); i++ ) 
            delete[] rend_path_strs[i];
    }

    delete[] lend_path_strs;
    delete[] rend_path_strs;


    lend_gsa = NULL;
    rend_gsa = NULL;
    lend_path_lens = NULL;
    rend_path_lens = NULL;
    lend_path_strs = NULL;
    rend_path_strs = NULL;

    gsa_ids.clear();

    lend_path_reads.clear();
    rend_path_reads.clear();

    lend_read_paths.clear();
    rend_read_paths.clear();
    
}
