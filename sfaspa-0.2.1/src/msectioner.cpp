#include "msectioner.h"

MergeSectioner::MergeSectioner( std::string *query, 
                                std::string *sbjct,
                                const KmerArray *qkmers,
                                const KmerArray *skmers,
                                const KmerPosMap *qp,
                                const KmerPosMap *sp )
{
    double t0 = mytime();
    
    Sectioner::init(query, sbjct);

    kmer_size = Param::merge_filter_kmer;
    extra     = Param::merge_expand_nbase;
    verbose   = Param::verbose;
    // aln_flag  = Param::merge_base_align;

    // merge_score = Param::merge_score;

    // align_score = 0.0;
    // align_success = false;



    // KmerArray query_kmers = biostr::getKmers(*query, kmer_size);
    // KmerArray sbjct_kmers = biostr::getKmers(*sbjct, kmer_size);

    query_kmers = qkmers;
    sbjct_kmers = skmers;

    qposs = qp;
    sposs = sp;
    /* beg = max, end = min */
    section = Section( query_kmers->size()-1+kmer_size, 0, sbjct_kmers->size()-1+kmer_size, 0, 0 );
    
    if ( verbose ) std::cout << "section init:" << mytime() -t0 << "\n";
}

bool 
MergeSectioner::find()
{
    double t0 = mytime();
    makePairMap();
    locate();
    //bool status = expand();
    expand();
    if ( verbose ) std::cout << "section find:" << mytime() -t0 << "\n";
    //return status;
    return true;
}

void MergeSectioner::makePairMap()
{
    double t0 = mytime();
    // std::tr1::unordered_map<KmerId, std::vector<int> > qposs, sposs;
    // for ( size_t i = 0; i < query_kmers.size(); i++ ) 
    //     qposs[query_kmers[i]].push_back(i);
    // for ( size_t i = 0; i < sbjct_kmers.size(); i++ ) 
    //     sposs[sbjct_kmers[i]].push_back(i);


    //std::tr1::unordered_map<KmerId, std::vector<int> >::iterator it;
    KmerPosMap::const_iterator it, jt;
    for ( it = sposs->begin(); it != sposs->end(); ++it ) {
        std::vector<size_t> sp = it->second;
        jt = qposs->find(it->first);
        if ( jt == qposs->end() ) continue;
        //if ( qposs->find(it->first) == qposs->end() ) continue;
        //std::vector<size_t> qp = (*qposs)[it->first];
        std::vector<size_t> qp = jt->second;
        for ( size_t i = 0; i < sp.size(); i++ ) 
            for ( size_t j = 0; j < qp.size(); j++ )
                poss_pair.insert(IntPair(qp[j], sp[i]));
    }
    if ( verbose ) std::cout << "make pair map:" << mytime() -t0 << "\n";
}

void MergeSectioner::locate()
{
    double t0 = mytime();

    PosPairMap::iterator qt;
    for ( qt = poss_pair.begin(); qt != poss_pair.end(); ++qt ) {
        //int diff = (qt->first-section.qbeg) - (qt->second-section.sbeg);

        if ( qt->first  < section.qbeg ) section.qbeg = qt->first;
        if ( qt->first  > section.qend ) section.qend = qt->first + kmer_size - 1 ;
        if ( qt->second < section.sbeg ) section.sbeg = qt->second;
        if ( qt->second > section.send ) section.send = qt->second + kmer_size - 1;
        section.count++;
    }
    int diff = (section.qend-section.qbeg) - (section.send-section.sbeg);
    if ( verbose ) 
        printf("region:\tqbeg:%d sbeg:%d qend:%d send:%d\tmatch:%d\tdiff:%d\n", section.qbeg,section.sbeg, section.qend,section.send, section.count, diff);

    // if ( section.count < mink ) return false;
    
    // return true;

    if ( verbose ) std::cout << "locate:" << mytime() -t0 << "\n";
}

// merge type
//bool 
void MergeSectioner::expand()
{
    //double t0 = mytime();

    // int qnkmers = (int)query_kmers.size();
    // int snkmers = (int)sbjct_kmers.size();

    int qsize = (int)query->size();
    int ssize = (int)sbjct->size();

    // section.sbeg ------------------------------ section.send (sbjct)
    //       section.qbeg ------------------ section.qend (query)

    // start expansion
    section.sbeg -= section.qbeg;
    section.qbeg = 0;
    if ( section.sbeg < 0 ) section.sbeg = 0;

    // end expansion
    // int soff = snkmers - (section.send+1);
    // int qoff = qnkmers - (section.qend+1);
    int soff = ssize - (section.send+1);
    int qoff = qsize - (section.qend+1);
    // section.qend = qnkmers-1;
    // section.send = (qoff<=soff) ? section.send+qoff : snkmers-1;
    section.qend = qsize-1;
    section.send = (qoff<=soff) ? section.send+qoff : ssize-1;

    if ( verbose ) 
        printf("expanded:\tqbeg:%d sbeg:%d qend:%d send:%d\n", section.qbeg,section.sbeg, section.qend,section.send);

    // if ( section.send-section.sbeg < section.qend-section.qbeg ) 
    //     return false;
    // return true;
}

//bool 
void MergeSectioner::pad()
{
    double t0 = mytime();
    //int qsize = (int)query->size();
    int ssize = (int)sbjct->size();

    // pad more
    int soff = ssize - (section.send+1);
    section.sbeg = ( section.sbeg >= extra ) ? section.sbeg-extra : 0;
    //section.send = ( soff >= extra ) ? section.send+extra : snkmers-1;
    section.send = ( soff >= extra ) ? section.send+extra : ssize-1;
 
    if ( verbose ) 
        printf("padded:\tqbeg:%d sbeg:%d qend:%d send:%d\n", section.qbeg,section.sbeg, section.qend,section.send);

    if ( section.send-section.sbeg < section.qend-section.qbeg ) {
        //section.sbeg = 0; section.send = snkmers-1;
        section.sbeg = 0; section.send = ssize-1;
        if ( verbose ) std::cout << "Can't handle. Use full range\n";
        // return false;
    }
    
    if ( verbose ) std::cout << "padded:" << mytime() -t0 << "\n";

    // return true;
}

