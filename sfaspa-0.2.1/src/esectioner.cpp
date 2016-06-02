#include "esectioner.h"

LatchSectioner::LatchSectioner( std::string *query, 
                                std::string *sbjct,
                                const KmerArray *qkmers,
                                const KmerArray *skmers,
                                const KmerPosMap *qp,
                                const KmerPosMap *sp)
{
    double t0 = mytime();

    Sectioner::init(query, sbjct);

    off_diag       = Param::extend_off_nbase;
    mink           = 1;
    overlap_length = Param::extend_length;

    anchor_kmer    = Param::extend_anchor_kmer;
    filter_kmer    = Param::extend_filter_kmer;
    filter_score   = Param::extend_filter_score;

    verbose        = Param::verbose;

    query_anchor_kmers = qkmers;
    sbjct_anchor_kmers = skmers;

    qposs = qp;
    sposs = sp;

    if ( verbose ) std::cout << "section init:" << mytime()-t0 << " sec\n";
}

bool LatchSectioner::find()
{
    return locate();
}

bool LatchSectioner::locate()
{
    double t0 = mytime();
    bool found = false;

    // region of interest
    // no. of k-mers to check
    int ncheck = overlap_length - anchor_kmer + 1;

    KmerArray check_kmers = direction == LEFT ? 
        KmerArray( sbjct_anchor_kmers->begin(),  sbjct_anchor_kmers->begin()  + ncheck ) :
        KmerArray( sbjct_anchor_kmers->rbegin(), sbjct_anchor_kmers->rbegin() + ncheck) ;
    
    // if ( Param::verbose ) {
    //     for ( size_t i  =0; i < check_kmers.size(); i++ )
    //         std::cout << alpha::IntegerToAminoAcid(check_kmers[i], anchor_kmer) << " ";
    //     std::cout << "\n";
    // }
    KmerArray::iterator it;
    KmerPosMap::const_iterator qt, st;
    for ( it = check_kmers.begin(); it != check_kmers.end(); ++it ) {
        qt = qposs->find(*it);
        if ( qt == qposs->end() ) continue;
        //if ( qposs->find(*it) == qposs->end() ) continue;

        std::vector<size_t>::const_iterator jt, kt;
        st = sposs->find(*it);        
        //for ( jt = (*sposs)[*it].begin(); jt != (*sposs)[*it].end(); ++jt ) {
        for ( jt = st->second.begin(); jt != st->second.end(); ++jt ) {
            int spos = *jt;
            //for ( kt = (*qposs)[*it].begin(); kt != (*qposs)[*it].end(); ++kt ) {
            for ( kt = qt->second.begin(); kt != qt->second.end(); ++kt ) {
                int qpos = *kt;
                if ( verbose ) std::cout << "spos:" << spos << "\tqpos:" << qpos << "\n";

                double tic = mytime();
                Section cs;
                found = makeRange(cs, spos, qpos);
                if ( verbose ) std::cout << "make range:" << mytime()-tic << "\n";
                if ( !found ) continue;
                
                tic = mytime();
                found = passKmerFilter(cs);
                if ( verbose ) std::cout << "kmer filter:" << mytime()-tic << "\n";
                if ( !found ) continue;
                
                if ( cs.count > section.count ) {
                    section = cs;
                    if ( verbose ) printf("New max region\n");
                }
            } 
        }
    }

    if ( verbose ) std::cout << "locate range:" << mytime()-t0 << "\n";
    return found;
}

bool LatchSectioner::makeRange( Section &sect, int spos, int qpos )
{
    int ss = spos;
    int se = spos + anchor_kmer - 1;    
    int qs = qpos;
    int qe = qpos + anchor_kmer - 1;

    if ( verbose ) printf("ss:%d se:%d qs:%d qe:%d\n", ss, se, qs, qe );
    int qend = (int)query->size()-1;
    int send = (int)sbjct->size()-1;
    
    int qoff = qend-qe;
    int soff = send-se;
    if ( verbose ) printf("soff:%d qoff:%d\n", soff, qoff);

    //=================================================
    //                  ss   se
    //             sbeg--ooooo-------------send (sbjct)
    //                   |||||
    // qbeg--------------ooooo--qend (query)    
    //                  qs   qe
    //=================================================
    if ( direction == LEFT ) {
        if ( qs <= ss ) {
            if ( verbose ) std::cout << "query inclusion\n";            
            return false;
        }
        if ( qoff >= soff ) {
            if ( verbose ) std::cout << "sbjct inclusion\n";
            return false;
        } 

        sect.sbeg = 0; 
        sect.qbeg = qs-ss;
        sect.qend = qend; 
        sect.send = se + qoff;
    }

    //=================================================
    //                  ss   se
    // sbeg--------------ooooo--send (sbjct)    
    //                   |||||
    //             qbeg--ooooo-------------qend (query)
    //                  qs   qe
    //=================================================
    else {
        if ( ss <= qs ) {
            if ( verbose ) std::cout << "sbjct inclusion\n";
            return false;
        }
        if ( soff >= qoff ) {
            if ( verbose ) std::cout << "query inclusion\n";
            return false;
        } 

        sect.qbeg = 0; 
        sect.sbeg = ss - qs;
        sect.send = send;
        sect.qend = qe + soff;
    }

    int diff = abs( (sect.qend-sect.qbeg) - (sect.send-sect.sbeg) );
    if ( diff > off_diag ) {
        if ( Param::verbose ) std::cout << "Off diagonal:" << diff << "\n";
        return false;
    }

    if ( verbose ) 
        printf("region:\tqs:%d ss:%d qe:%d se:%d\n", sect.qbeg,sect.sbeg,sect.qend,sect.send);

    return true;
}


//====================================================================
// 
//====================================================================
bool LatchSectioner::passKmerFilter( Section &s )
{
    int len_ext = (s.qend-s.qbeg < s.send-s.sbeg) ? 
        (s.qend-s.qbeg+1) : (s.send-s.sbeg+1);
    assert(len_ext>0);

    if ( verbose ) printf("filter-kmer:%d\tfilter-score:%.2f\n", 
                          filter_kmer, filter_score );

    int min_ext = filter::minSameKmerCount( len_ext, filter_kmer, 1-filter_score ) ;
    if ( min_ext < 1 ) return false;

    std::string ssub = sbjct->substr(s.sbeg, s.send-s.sbeg+1);
    std::string qsub = query->substr(s.qbeg, s.qend-s.qbeg+1);

    if ( Param::verbose ) {
        std::cout << "ssub:" << ssub << "\n";
        std::cout << "qsub:" << qsub << "\n";
    }

    if ( (int)ssub.size() < filter_kmer || (int)qsub.size() < filter_kmer ) return false;

    KmerArray query_filter_kmers = biostr::getKmers(qsub, filter_kmer);
    KmerArray sbjct_filter_kmers = biostr::getKmers(ssub, filter_kmer);

    std::set<KmerId> query_set;
    for ( size_t a = 0; a < query_filter_kmers.size(); a++ )
        query_set.insert( query_filter_kmers[a] );

    int num_ext = 0;
    for ( size_t b = 0; b < sbjct_filter_kmers.size(); b++ )
        if ( query_set.find(sbjct_filter_kmers[b]) != query_set.end() )
            num_ext++;
    
    if ( verbose ) printf("extended region:\tqs:%d ss:%d qe:%d se:%d\n", s.qbeg,s.sbeg,s.qend,s.send);
    if ( verbose ) std::cout << "extended length:" << len_ext << "\tmin k-mer match:" << min_ext << "\t#match:" << num_ext << "\n";

    if ( num_ext < min_ext ) return false;

    s.count = num_ext;

    return true;        
}

