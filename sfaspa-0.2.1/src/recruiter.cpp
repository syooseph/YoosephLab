#include "recruiter.h"

Recruiter::Recruiter()
{

}

Recruiter::~Recruiter()
{
    clear();
}

void Recruiter::clear()
{
    kmer_paths.clear();
    path_kposs.clear();
}

void Recruiter::run()
{
    derived = DERIVED_RECRUIT;
    recruit();
    write();

    status = true;
}


void Recruiter::recruit()
{
    updatePathMembership();

    makeIndex();

    ReadIdList rlist;
    if ( Param::verbose >= 1 ) 
        std::cout << "# total reads:" << nreads << "\n";
    assert( preads != NULL );
    for ( int k = 0; k < nreads; k++ ) 
        if ( preads[k] == NOT_PATH && preads[k] != BAD_READ ) rlist.push_back((ReadId)k);

    ReadIdArray orphans = ReadIdArray( rlist.begin(), rlist.end() );
    if ( Param::verbose >= 1 ) 
        std::cout << "# Unrecruited:" << orphans.size() << "\n";

    //Param::ncpus == 1 ? __recruitSP(orphans) : __recruitMP(orphans);
    __recruit(orphans);

    int assem = MetaPath::countAssembledReads();
    double ratio = (double)assem/nreads;
    printf("assembled read ratio: assembled/total-reads = %.2f%% (%d/%d)\n", ratio*100, assem, nreads);
}

// void Recruiter::__recruitSP( ReadIdArray &orphans )
// {
//     //double t0 = mytime();
//     size_t i, n = orphans.size();
//     progress = Progress( 1.0, 0, n, mytime() );

//     for ( i = 0; i < n; i++ ) {
//         ReadPlacement place;
//         ReadId rid = orphans[i];
//         PathId pid = findPath( rid, place );
//         if ( preads[rid] == BAD_READ ) continue;
//         if ( pid != NOT_PATH ) {            
//             size_t nid = Path2Ids[pid];
//             Alns[nid].addPlacement(place);
//             preads[rid] = pid;
//         }

//         double prev = progress.ratio;
//         progress.count++;
//         progress.showProgress();
    
//         if ( progress.ratio > prev && ( Param::summary_flag || Param::verbose ) ) {
//             recruit_log.print();//printSummary(t0);
//             anchor_log.print();
//         }
//     }

//     if ( !Param::summary_flag && !Param::verbose )
//         recruit_log.print();//printSummary(t0);
//         //printSummary(t0);
// }

void Recruiter::__recruit( ReadIdArray &orphans )
{
    size_t i, n = orphans.size();
    progress = Progress( 1.0, 0, n, mytime() );

    bool par_run = Param::ncpus > 1 ? true : false;
#pragma omp parallel for if (par_run) schedule(dynamic, 1) private(i) num_threads(Param::ncpus)    
    for ( i = 0; i < n; i++ ) {
        ReadPlacement place;
        ReadId rid = orphans[i];
        if ( preads[rid] == BAD_READ ) continue;

        PathId pid = findPath( rid, place );
        if ( pid != NOT_PATH ) {            
            if ( Param::verbose > 1 ) std::cout << "Success\trid:" << rid << "\tpid:" << pid << "\n";
            size_t nid = Path2Ids[pid];
#pragma omp critical
{
            Alns[nid].addPlacement(place);
            Alns[nid].setDirty(true);
}
            preads[rid] = pid;
        }

        double prev = progress.ratio;
#pragma omp critical /* single might be better */
{
        progress.count++;
        progress.showProgress();
    
        if ( progress.ratio > prev && Param::verbose >= 1 ) 
            recruit_log.print();
} 
    }
    if ( ! Param::verbose ) 
        recruit_log.print();
}

void Recruiter::makeIndex()
{
    double t0 = mytime();
    for ( size_t i = 0; i < npaths; i++ ) {
        if ( bad_paths[i] ) continue;
        assert( Id2Paths.find(i) != Id2Paths.end() );
        PathId pid = Id2Paths[i];
        std::string str = Alns[i].getSequence();
        if ( Param::verbose > 1 ) 
            std::cout << i << "\tpath:" << pid<< "\tbad:" << bad_paths[i] << "\t" << str << "\n";
        assert(str.size());
        KmerArray kmers = biostr::getKmers(str, Param::recruit_filter_kmer);
        assert(kmers.size());

        KmerPosMap kposs;
        for ( size_t j = 0; j < kmers.size(); j++ )
            kposs[kmers[j]].push_back(j);

        
        //path_kposs.insert( std::pair<PathId, KmerArray>( pid, kmers ) );
        //path_kposs.add(pid, kmers);
        path_kposs.add(pid, kposs);

        kmer_paths.add(kmers, pid);
    }
    if ( Param::verbose >= 1 ) 
        std::cout << "Kmer to path mapping initialized:" << mytime()-t0 << " sec\n";
}

PathId Recruiter::findPath( ReadId rid,
                         ReadPlacement &place )
{
    std::string read = seqs[rid];
    assert(read.size() > 0);
    if ( Param::verbose > 1 ) {
#pragma omp critical
        std::cout << "Rid:" << rid << "\t" << read << "\n";
    }
    
    KmerArray query_kmers = biostr::getKmers( read, Param::recruit_filter_kmer );

    if ( query_kmers.size() == 0 ) return NOT_PATH;

    int min_filter = filter::minSameKmerCount( read.size(), Param::recruit_filter_kmer, 1-Param::recruit_filter_score );
    if ( min_filter < Param::recruit_min_filter ) 
        min_filter = Param::recruit_min_filter;

    if (Param::verbose > 1) {
#pragma omp critical
        std::cout << "Min k-mer:" << min_filter << "\n";
    }

    std::multimap<size_t, PathId> freq_map;
    findCandidatePaths(read, query_kmers, freq_map, min_filter);
    if ( freq_map.size() == 0 ) return NOT_PATH;

    return findBestPathByAnchoring(rid, read, freq_map, place, query_kmers, min_filter);
    // return Param::recruit_by_align ?
    //     findBestPathByAlignment(rid, read, freq_map, place, query_kmers, min_filter) :
    //     findBestPathByAnchoring(rid, read, freq_map, place, query_kmers, min_filter);
}

void Recruiter::findCandidatePaths( std::string &read,
                                    KmerArray &query_kmers,
                                    std::multimap<size_t, PathId> &freq_map,
                                    int min_filter
                                    )
{
    double t0 = mytime();        


    double tic = mytime();
    PathFreqMap path_freq;
    kmer_paths.findPaths( path_freq, query_kmers );
    recruit_log.et_can_all += (mytime()-tic);
    // std::multimap<size_t, PathId> path_map = util::sortByValue<PathId, size_t>(path_freq);
    // std::multimap<size_t, PathId>::reverse_iterator it;
    // for ( it = path_map.rbegin(); it != path_map.rend(); ++it) {
    //     if ( (int)it->first < min_filter ) break;
    //     freq_map.insert(std::pair<size_t, PathId>(it->first, it->second));
    // }
    
    tic = mytime();
    //for ( PathFreqMap::iterator it = path_freq.begin(); it != path_freq.end(); ++it ) {
    for ( auto entry : path_freq ) {
        //if ( (int)it->second < min_filter ) continue;
        if ( (int)entry.second < min_filter ) continue;
        //freq_map.insert(std::pair<size_t, PathId>(it->second, it->first));
        freq_map.insert(std::pair<size_t, PathId>(entry.second, entry.first));
    }

    if ( Param::verbose > 1) {
#pragma omp critical
        std::cout << "#Paths:" << freq_map.size() 
                  << "\tScan paths time:" << mytime() -t0 << "\n";
    }
    recruit_log.ct_can += freq_map.size();
    recruit_log.et_can_good += (mytime()-tic);
    recruit_log.et_can += (mytime()-t0);
}

void Recruiter::makePosPairMap( std::multimap<int,int> &poss_pair,
                                KmerArray &query_kmers,
                                KmerArray &sbjct_kmers )
{
    std::tr1::unordered_map<KmerId, std::vector<int> > qposs, sposs;
    for ( size_t i = 0; i < query_kmers.size(); i++ ) 
        qposs[query_kmers[i]].push_back(i);
    for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) 
        sposs[sbjct_kmers[i]].push_back(i);


    std::tr1::unordered_map<KmerId, std::vector<int> >::iterator it;
    for ( it = sposs.begin(); it != sposs.end(); ++it ) {
        std::vector<int> sp = it->second;
        if ( qposs.find(it->first) == qposs.end() ) continue;
        std::vector<int> qp = qposs[it->first];
        for ( size_t i = 0; i < sp.size(); i++ ) 
            for ( size_t j = 0; j < qp.size(); j++ )
                poss_pair.insert(IntPair(qp[j], sp[i]));
    }
}

bool Recruiter::findRegion( std::multimap<int,int> &poss_pair,
                            std::multimap<int,int>::iterator &pt,
                            IntPair &rb,
                            IntPair &sb,
                            int &ct,
                            KmerArray &query_kmers )
{
    rb.first = rb.second = pt->first;  // query
    sb.first = sb.second = pt->second; // sbjct
    ct = 0;
    
    if ( Param::verbose > 1 ) {
        std::cout << "Kmer:" << alpha::IntegerToAminoAcid(query_kmers[rb.first], Param::extend_filter_kmer) << "\n";
        printf("init\ts1:%d s2:%d e1:%d e2:%d\n", rb.first,sb.first,rb.second,sb.second);
    }

    std::multimap<int,int>::iterator qt;
    for ( qt = pt; qt != poss_pair.end(); ++qt ) {
        if ( Param::verbose > 1 ) {
            std::cout << "Kmer:" << alpha::IntegerToAminoAcid(query_kmers[qt->first], Param::extend_filter_kmer) << "\n";
            printf("loop\ts1:%d s2:%d e1:%d e2:%d\n", rb.first,sb.first, qt->first, qt->second);
        }
        if ( qt->first  < rb.first  || qt->first  < rb.second ) continue;
        if ( qt->second < sb.first  || qt->second < sb.second ) continue;
        int diff = (qt->first-rb.first) - (qt->second-sb.first);
        if ( abs(diff) > 1 ) continue;
        rb.second = qt->first;
        sb.second = qt->second;
        ct++;
    }
    
    if ( Param::verbose > 1 ) printf("out\trb.first:%d s2:%d e1:%d e2:%d\tmatch:%d\n", rb.first,sb.first, rb.second,sb.second, ct);
    int diff = (rb.second-rb.first) - (sb.second-sb.first);
    if ( Param::verbose > 1 ) std::cout << "diff:" << diff << "\n";
    if ( abs(diff) > 1 ) return false;
    
    return true;
}

void Recruiter::extendRegion( IntPair &rb,
                              IntPair &sb,
                              int qnkmers,
                              int snkmers )
{
    if ( Param::verbose > 1 ) 
    printf("rs:%d,re:%d,ss:%d,se:%d\n", rb.first,rb.second,sb.first,sb.second);

    if ( rb.first > 0 ) {
        if ( sb.first > rb.first ) { sb.first -= rb.first; rb.first = 0; }
        else { rb.first -= sb.first; sb.first = 0; }
    }
    if ( qnkmers-rb.second > 1 ) {
        if ( snkmers-sb.second > qnkmers-rb.second ) {
            sb.second += (qnkmers-rb.second-1);
            rb.second = qnkmers-1;
        }
        else {
            rb.second += (snkmers-sb.second-1);
            sb.second = snkmers-1;
        }
    }
    if ( Param::verbose > 1 ) 
    printf("rs:%d,re:%d,ss:%d,se:%d\n", rb.first,rb.second,sb.first,sb.second);
}

bool Recruiter::passKmerFilter( IntPair &rb,
                                IntPair &sb,
                                KmerArray &query_kmers,
                                KmerArray &sbjct_kmers )
{
    int len_ext = (rb.second-rb.first < sb.second-sb.first) ? 
        (rb.second-rb.first+Param::recruit_filter_kmer) :
        (sb.second-sb.first+Param::recruit_filter_kmer) ;
    
    int min_ext = filter::minSameKmerCount( len_ext, Param::recruit_filter_kmer, 1-Param::recruit_filter_score );
    int num_ext = 0;
    
    std::set<KmerId> query_set;
    for ( int a = rb.first; a <= rb.second; a++ ) 
        query_set.insert( query_kmers[a] );
    for ( int b = sb.first; b <= sb.second; b++ )
        if ( query_set.find(sbjct_kmers[b]) != query_set.end() )
            num_ext++;
    
    if ( Param::verbose > 1 ) printf("ext\ts1:%d s2:%d e1:%d e2:%d\n", rb.first,sb.first, rb.second, sb.second);
    if ( Param::verbose > 1 ) std::cout << "recruit length:" << len_ext << "\tmin k-mer:" << min_ext << "\t#ext-match:" << num_ext << "\n";

    if ( num_ext < min_ext ) return false;
    return true;        
}

IntPair Recruiter::findSbjctRegion(std::string &read, 
                                   std::string &sbjct,
                                   KmerArray &query_kmers,
                                   KmerArray &sbjct_kmers)
{
    int min = filter::minSameKmerCount( read.size(), Param::recruit_filter_kmer, 1-Param::recruit_filter_score );
    if ( min < 1 ) min = 1;

    // KmerArray query_kmers = biostr::getKmers( read , Param::recruit_filter_kmer );
    // KmerArray sbjct_kmers = biostr::getKmers( sbjct, Param::recruit_filter_kmer );

    std::multimap<int,int> poss_pair;
    makePosPairMap(poss_pair, query_kmers, sbjct_kmers );

    IntPair rb, sb, rbmax, sbmax(0, sbjct.size()-1);
    int ct, max = 0;
    std::multimap<int,int>::iterator pt,qt;
    for ( pt = poss_pair.begin(); pt != poss_pair.end(); ++pt ) {

        if ( ! findRegion(poss_pair,pt,rb,sb,ct,query_kmers) ) continue;

        if ( Param::verbose > 1 ) {
            std::cout << "#kmer-match:" << ct << "\n";
        }

        if ( ct < min ) continue;
        if ( ct < max ) continue;                            
        if ( ct > max ) {
            extendRegion( rb,sb, query_kmers.size(), sbjct_kmers.size() );
            
            if ( passKmerFilter( rb, sb, query_kmers, sbjct_kmers ) ) {
                sb.second += (Param::recruit_filter_kmer-1);
                sbmax = sb;
                max = ct;
            } else continue;
        }
    }
    return sbmax;
}


PathId Recruiter::findBestPathByAnchoring( ReadId rid,
                                           std::string &read,
                                           std::multimap<size_t, PathId> &freq_map,
                                           ReadPlacement &place,
                                           KmerArray &query_kmers,
                                           int min_filter)
{
    if ( freq_map.size() == 0 ) return NOT_PATH;

    size_t niter = 0;
    std::multimap<size_t, PathId>::reverse_iterator it;
    for ( it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
        double t0 = mytime();

#pragma omp atomic
        recruit_log.ct_cmp++;

        PathId pid = it->second;
        assert( Path2Ids.find(pid) != Path2Ids.end() );
        size_t i = Path2Ids[pid];
        std::string sbjct = Alns[i].getSequence();
        assert(sbjct.size()>0);

        recruit_log.et_sbjct_find += (mytime()-t0);

        if ( sbjct.size() < read.size() ) return NOT_PATH;
        
        if ( Param::verbose > 1 ) 
#pragma omp critical
            std::cout << "Iteration:" << ++niter 
                      << "\tPath:" << pid
                      << "\tLength:" << sbjct.size() 
                      << "\t# Shared kmers:" << it->first << "\n"
                      << sbjct << "\n"
                      << mytime()-t0 << " sec\n";
        
        //assert( path_kposs.find(path) != path_kposs.end() );
        //assert( path_kposs.has(path) );
        //assert( path_kposs[path].size() > 0 );

        double tic = mytime();
        //KmerArray sbjct_kmers = path_kposs.getPointer(pid);
        const KmerPosMap* sbjct_kposs = path_kposs.getPointer(pid);
        // std::string sbj =  biostr::getSequenceString( &sbjct_kmers[0], sbjct_kmers.size(), Param::recruit_filter_kmer );
        // if ( sbjct != sbj ) {
        //     std::cout << "Sbjct sequence match error\n"
        //               << "sbjct:" << pid << "\n"
        //               << "sbjct1:" << sbjct << "\n"
        //               << "sbjct2:" << sbj << "\n";
        // }
        // assert( sbjct == sbj );

        assert( sbjct_kposs->size() > 0 );
        
        recruit_log.et_sbjct_get += (mytime()-tic);

        int k = Param::recruit_filter_kmer;
        //bool v = Param::verbose;
        bool v = Param::verbose > 1 ? 1 : 0;
        tic = mytime();
        //Anchor anchor( &sbjct, &read, &(path_kposs[path]), &query_kmers, k, v );
        //Anchor anchor( rid, path, &sbjct, &read, &sbjct_kmers, &query_kmers, k, v );
        //Anchor anchor( rid, path, sbjct, read, sbjct_kmers, query_kmers, k, v );
        //Anchor *anchor = new Anchor( rid, pid, sbjct, read, sbjct_kmers, query_kmers, k, v );
        //Anchor anchor( rid, pid, &sbjct, &read, &sbjct_kmers, &query_kmers, min_filter, k, v );
        Anchor anchor( rid, pid, &sbjct, &read, sbjct_kposs, &query_kmers, min_filter, k, v );
        bool found = anchor.find();
        
        AnchorLog alog = anchor.getLog();        
        anchor_log.add(alog);
        
        // KmerPossMap sposs;
        // for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) 
        //     sposs[sbjct_kmers[i]].push_back(i);

        // KmerPossMap qposs;
        // for ( size_t i = 0; i < query_kmers.size(); i++ ) 
        //     qposs[query_kmers[i]].push_back(i);

        // MergeSectioner sect(&read, &sbjct, &query_kmers, &sbjct_kmers, &qposs, &sposs);
        // Section s = sect.getSection();

        if ( Param::verbose > 1 ) {
#pragma omp critical
            std::cout << "Anchor search:" << mytime()-t0 << " sec\n";
            std::cout << "Result:" << found << "\n";
        }
        recruit_log.et_anc += (mytime()-tic);
        
        // if ( s.send-s.sbeg != s.qend-s.qbeg) continue;

        // std::string sbjct_sub = sbjct.substr(s.sbeg, s.send-s.sbeg+1);
        // std::string query_sub = read.substr(s.qbeg, s.qend-s.qbeg+1);
    
        
        //if ( !anchor.getStatus() ) {
        if ( !found ) {
            if ( Param::verbose > 1 ) {
#pragma omp critical
                std::cout << "Anchoring fail\n";
            }
            //delete anchor;
            continue;
        }

        
        // size_t length = anchor.getLength();
        // if ( (double)length/read.size() < Param::recruit_ratio ) {
        //     //delete anchor;
        //     continue;
        // }
        // // if ( (double)query_sub.size()/read.size() < Param::recruit_ratio ) 
        // //     continue;

        double score = anchor.getScore();
        if ( Param::verbose > 1 ) std::cout << "Score:" << score << "\n";
        if ( score < Param::recruit_score ) {
            //delete anchor;
            continue;
        }

        // double score = scoring::positiveRate( sbjct_sub, query_sub, BLOSUM62 );
        // if ( score < Param::recruit_score ) {
        //     continue;
        // }
        
        // IntPair sposs = anchor.getSbjctPos();
        // IntPair qposs = anchor.getQueryPos();

        Section s = anchor.getSection();

        if ( Param::verbose > 1 ) {
#pragma omp critical
            //std::cout << "region:" << qposs.first << " " << qposs.second << "\t" << mytime()-t0 << " sec\n";
            std::cout << "region:" << s.qbeg << " " << s.qend << "\t" << mytime()-t0 << " sec\n";
            }
         place.rid = rid;
         // place.read_pos = qposs.first; 
         // place.path_pos = sposs.first; 
         //place.length = length;
         assert(s.qbeg == 0 );
         place.read_pos = 0; 
         place.path_pos = s.sbeg; 
         place.length = read.size();
             
#pragma omp atomic
        recruit_log.ct_good++;

        //delete anchor;
        return pid;
    }

    return NOT_PATH;
}


// PathId Recruiter::findBestPathByAlignment( ReadId rid,
//                                            std::string &read,
//                                            std::multimap<size_t, PathId> &freq_map,
//                                            ReadPlacement &place,
//                                            KmerArray &query_kmers,
//                                            int min_filter)
// {
//     if ( freq_map.size() == 0 ) return NOT_PATH;
    
//     size_t niter = 0;
//     std::multimap<size_t, PathId>::reverse_iterator it;
//     for ( it = freq_map.rbegin(); it != freq_map.rend(); ++it ) {
//         PathId path = it->second;
//         size_t i = Path2Ids[path];
//         std::string sbjct = Alns[i].getSequence();
//         assert(sbjct.size()>0);

//         if ( sbjct.size() < read.size() ) return NOT_PATH;

//         if ( Param::verbose ) 
//             std::cout << "Iteration:" << ++niter 
//                       << "\tPath:" << path
//                       << "\tLength:" << sbjct.size() 
//                       << "\t# Shared kmers:" << it->first << "\n"
//                       << sbjct << "\n";

//         double t0 = mytime();

//         //assert( path_kposs.find(path) != path_kposs.end() );

//         KmerArray sbjct_kmers = path_kposs.getPointer(path);

//         IntPair poss = IntPair(0, sbjct.size()-1);
//         if ( sbjct.size() >= 2*read.size() ) {
//             //poss = findSbjctRegion(read, sbjct, query_kmers, path_kposs [path]);
//             poss = findSbjctRegion(read, sbjct, query_kmers, sbjct_kmers);
//             assert(poss.first>=0);
//             assert(poss.first<(int)sbjct.size());
//             assert(poss.second<(int)sbjct.size());
//         }

//         recruit_log.et_pos += (mytime()-t0);

//         if ( Param::verbose ) 
//             std::cout << "region:" << poss.first << " " << poss.second << "\t" << mytime()-t0 << " sec\n";
        
//         int off = poss.first;
//         int len = poss.second-poss.first+1;
//         assert(off+len <= (int)sbjct.size());
//         std::string sbjct_sub = sbjct.substr(off,len);
//         if ( Param::verbose ) std::cout << "sub:" << sbjct_sub << "\n";
//         if ( sbjct_sub.size() >= read.size() ) {
//             t0 = mytime();
//             size_t found = sbjct_sub.find(read);
//             recruit_log.et_sub += (mytime()-t0);

// #pragma omp atomic 
//             recruit_log.ct_sub++;

//             if ( found != std::string::npos ) {
//                 place.rid = rid; place.read_pos = 0; place.path_pos = found+off;
//                 place.length = read.size();
// #pragma omp atomic
//                 recruit_log.sub_good++;
//                 return path;
//             } 
//         } 
        
//         t0 = mytime();
//         SemiGlobalAlign aln(sbjct_sub, read, ANCHOR_CENTER, Param::gap_ext, Param::gap_open);
//         recruit_log.et_aln += (mytime()-t0);
// #pragma omp atomic
//         recruit_log.ct_aln++;
//         AlignSummary summary = aln.getSummary();
        
//         if ( Param::verbose ) aln.printAlignment(std::cout);
//         if ( Param::verbose ) summary.print(std::cout);
        
//         double score = !Param::identity_flag ?
//             summary.posrate : 
//             (double)summary.match/summary.length;
//         if ( score < Param::recruit_score ) continue;
//         if ( (double)summary.length/read.size() < Param::recruit_ratio ) continue;
//         if ( summary.lgap.first && summary.lgap.second ) continue;
        
//         place.rid = rid;
//         if ( summary.lgap.first > 0 ) {
//             place.read_pos = summary.lgap.first; place.path_pos = off; place.length = summary.length;
//         } else {
//             place.read_pos = 0; place.path_pos = off+summary.lgap.second; place.length = summary.length;
//             AlignPosList::iterator it;
//             for (  it = summary.ilist.begin(); it != summary.ilist.end(); ++it )
//                 place.ilist.push_back( Mismatch(rid, it->qry_pos, off+it->ref_pos ) );
//             for (  it = summary.dlist.begin(); it != summary.dlist.end(); ++it )
//                 place.dlist.push_back( Mismatch(rid, it->qry_pos, off+it->ref_pos ) );
//         }
// #pragma omp atomic
//         recruit_log.aln_good++;
//         return path;
//     }
    
//     return NOT_PATH;
// }

// void Recruiter::printSummary( double t0)
// {
//     printf("# Candidates: %zu\n", ct_can);
//     if ( ! Param::recruit_by_align ) {
//         printf("# Comparisons: %zu\n", ct_cmp );
//         printf("# Recruited: %zu\n", ct_good);
//         printf("# Elapsed:%.2f sec (search:%.2f [paths:%.2f, trim:%.2f], anchoring:%.2f, sbjct-find:%.2f sbjct-get:%.2f)\n", mytime()-t0, et_can, et_can_all, et_can_good, et_anc, et_sbjct_find, et_sbjct_get);
//     } else {
//         printf("# Comparisons: %zu (substring:%zu, alignment:%zu)\n", 
//                ct_sub+ct_aln, ct_sub, ct_aln );
//         printf("# Recruited: %zu (substring:%zu, alignment:%zu)\n",
//                sub_good+aln_good, sub_good, aln_good);
        
//         printf("# Elapsed:%.2f (search:%.2f, substring:%.2f, range-search:%.2f, alignment:%.2f)\n",
//                mytime()-t0, et_can, et_sub, et_pos, et_aln);
//     }

// }
