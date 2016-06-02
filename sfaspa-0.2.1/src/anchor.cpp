#include "anchor.h"

Anchor::Anchor( std::string *s, 
                std::string *q, 
                const KmerPosMap *sk, 
                const KmerArray *qk,
                int m,
                int k, 
                bool v )
{
    sbjct      = s;
    query      = q;
    kmer       = k;
    verbose    = v;
    skposs     = sk;
    qkmers     = qk;
    min_filter = m;

    found      = false;
}

Anchor::Anchor( ReadId &r,
                PathId &p,
                std::string *s, 
                std::string *q, 
                const KmerPosMap *sk, 
                const KmerArray *qk,
                int m,
                int k, 
                bool v )
{
    query_rid  = r;
    sbjct_pid  = p;
    sbjct      = s;
    query      = q;
    kmer       = k;
    verbose    = v;
    skposs     = sk;
    qkmers     = qk;
    min_filter = m;

    found      = false;
}

Anchor::~Anchor()
{

}


//====================================================================
// Fine alignment region and calculate simiarity of the region by
// comparing base by base.
//====================================================================
bool Anchor::find()
{
    initPosMap();
    return findAnchor();
}


//====================================================================
// For each kmer in query, find kmer positions in sbjct sequence.
//====================================================================
void Anchor::initPosMap()
{
    double t0 = mytime();

    KmerPosMap::const_iterator it;
    for ( size_t qpos = 0; qpos < qkmers->size(); qpos++ ) {
        it = skposs->find( (*qkmers)[qpos] );
        if ( it == skposs->end() ) continue;

        for ( auto spos : it->second ) 
            poss_pair.insert( IntPair( (int)qpos, (int)spos ) );
    }
    log.et_init = mytime()-t0;
}

//====================================================================
// Find a good anchor region
//====================================================================
bool Anchor::findAnchor()
{
    double t0 = mytime();

    setbuf(stdout, NULL);
    found = false;

    for ( PosPairMap::iterator it = poss_pair.begin(); it != poss_pair.end(); ++it ) {
        log.ct_all++;
        
        //----------------------
        // Find an anchor region
        //----------------------
        double tic = mytime();
        bool good_range = findRange( it->first, it->second );
        log.et_range += (mytime()-tic);

        if ( !good_range ) continue;
        if ( verbose ) std::cout << "Found:" << found << "\n";
        
        log.ct_valid++;

        //------------------------------------
        // If already checked, skip this range
        //------------------------------------
        if ( exist() ) {
            log.ct_exist++;
            if ( verbose ) std::cout << "Already exist\n";
            continue;
        }

        //-------------------------------------------------------
        // Insert this one to checked list so as not to use later
        //-------------------------------------------------------
        tic = mytime();
        checked.insert( std::pair<IntPair, bool> ( IntPair(section.sbeg, section.send), true ) );
        log.et_map += (mytime()-tic);

        //----------------------------
        // Compute matching kmer count
        //----------------------------
        int match = 0;
        for ( auto kmer : *qkmers ) {
            if ( skposs->find(kmer) == skposs->end() ) continue;
            match++;
        }

        //-------------------------------------------
        // Should have at least min_filter kmer match
        //-------------------------------------------
        if ( match < min_filter ) continue;
        if ( verbose ) std::cout << "match:" << match << "\tmin-filter:" << min_filter << "\n";
        log.ct_pass++;

        length = section.qend - section.qbeg + 1;
        assert( length == query->size() );

        //------------------------------------------
        // Get a substring of sbjct in anchor region
        //------------------------------------------
        tic = mytime();
        std::string strim = sbjct->substr( section.sbeg, length );
        log.et_substr += (mytime()-tic);

        //---------------------------------------------------------------
        // Compute similarity score between query and the pivot substring
        // Compare base by base using BLOSUM62 matrix.
        //---------------------------------------------------------------
        tic = mytime();
        score = scoring::positiveRate( strim, *query, BLOSUM62 );
        log.et_score += (mytime()-tic);

        if ( verbose ) {
            std::cout << "strim:" << strim << "\n";
            std::cout << "score:" << score << "\n";
        }

        //--------------------------------
        // Break the loop if score is good
        //--------------------------------
        if ( score >= Param::recruit_score ) {
            found = true;
            break;
        }
    }

    log.et_find = mytime()-t0;
    return found;
}

//====================================================================
// Extract an anchor region.
// Starting from a single kmer point as start, extend until a path
// reaches to end. 
// When all the query sequence is region is covered, it is considered
// as good anchor region. 
// If sbjct sequence ends before query sequence, it is considered as
// a bad one.
//====================================================================
bool Anchor::findRange(int qpos, int spos )
{
    section.clear();

    //------------------------------------------
    // Initial start/end position of sbjct/query
    //------------------------------------------
    int ss = spos;
    int se = spos + kmer - 1;    
    int qs = qpos;
    int qe = qpos + kmer - 1;

    if ( verbose ) printf("ss:%d se:%d qs:%d qe:%d\n", ss, se, qs, qe );
    int qend = (int)query->size()-1;
    int send = (int)sbjct->size()-1;
    
    int qoff = qend-qe;
    int soff = send-se;
    if ( verbose ) printf("soff:%d qoff:%d\n", soff, qoff);

    //=================================================
    //         ss             se
    // sbeg-----oxooooooooooooo-------------send (sbjct)
    //          | |||||||||||||
    //    qbeg--oxooooooooooooo--qend (query)    
    //         qs             qe
    //=================================================
    
    if ( qs > ss ) return false;
    if ( qoff > soff ) return false;
    
    //-----------------------------------------------------------
    // Extend initial anchor region to cover entire query region.
    // Then, save it to range section.
    //-----------------------------------------------------------
    section.sbeg = ss - qs;
    section.qbeg = 0;
    section.qend = qend; 
    section.send = se + qoff;

    //-------------------------------------------------
    // Length of two sequence regions must be identical
    //-------------------------------------------------
    if ( section.send-section.sbeg != section.qend-section.qbeg ) return false;

    if ( verbose ) printf("ss:%d se:%d qs:%d qe:%d\n", section.sbeg, section.send, section.qbeg, section.qend);

    return true;
}

//====================================================================
// Check whether a current region was previously checked one.
//====================================================================
bool Anchor::exist()
{
    IntPair sbj(section.sbeg, section.send);
    return checked.find(sbj) != checked.end();
}

