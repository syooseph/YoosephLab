#include "path.h"

SpaPath::SpaPath()
{
    nkmer = nread = cins = cdel = 0;
    kmers = NULL;
    reads = NULL;
    inits = NULL;
    inss  = NULL;
    dels  = NULL;
    vprof = NULL;
    consensus = NULL;
}

SpaPath::SpaPath( KmerId *kids, ReadId *rids, unsigned nk, unsigned nr )
{
    build( kids, rids, nk, nr );
}

SpaPath::SpaPath( const char *seq, KmerId *kids, ReadId *rids, unsigned nk, unsigned nr )
{
    int len = strlen(seq);
    if ( len ) {
        consensus = new char[len+1];
        strncpy(consensus, seq, len+1);
    }
    build( kids, rids, nk, nr );
}

void SpaPath::build(KmerId *kids, ReadId *rids, unsigned nk, unsigned nr)
{
    nkmer = nk;
    nread = nr; 
    
    if ( nkmer > 0 ) kmers = new KmerId[nkmer];
    if ( nread > 0 ) reads = new ReadId[nread];
    if ( nread > 0 ) inits = new int[nread];

    for ( unsigned i = 0; i < nk; i++ ) kmers[i] = kids[i];
    for ( unsigned i = 0; i < nr; i++ ) reads[i] = rids[i];

    cins = cdel = 0;
    
    vprof = new ProfileVector();
}

SpaPath::SpaPath( const SpaPath &source )
{
    __copy(source);
}

SpaPath& SpaPath::operator= ( const SpaPath &source )
{
    __copy(source);
    return *this;
}

void SpaPath::__copy( const SpaPath &source )
{
    nkmer = source.nkmer;
    nread = source.nread;
    if ( nkmer > 0 ) kmers = new KmerId[nkmer];
    if ( nread > 0 ) reads = new ReadId[nread];
    if ( nread > 0 ) inits = new int[nread];
    for ( unsigned i = 0; i < nkmer; i++ ) kmers[i] = source.kmers[i];
    for ( unsigned i = 0; i < nread; i++ ) reads[i] = source.reads[i];
    for ( unsigned i = 0; i < nread; i++ ) inits[i] = source.inits[i];
    
    int nstr = strlen(source.consensus);
    if ( nstr > 0 ) {
        consensus = new char[nstr+1];
        strcpy(consensus, source.consensus);
    }

    cins = source.cins;
    cdel = source.cdel;
    if ( cins > 0 ) inss = new Mismatch[cins];
    if ( cdel > 0 ) dels = new Mismatch[cdel];
    for ( unsigned i = 0; i < cins; i++ ) inss[i] = source.inss[i];
    for ( unsigned i = 0; i < cdel; i++ ) dels[i] = source.dels[i];
    
    vprof = new ProfileVector(*source.vprof);
    //vprof = source.vprof;
}

SpaPath::~SpaPath()
{ 
    clear(); 
}

void SpaPath::clear()
{
    //if ( consensus != NULL && strlen(consensus) > 0 ) delete[] consensus;
    if ( strlen(consensus) > 0 ) delete[] consensus;
    if ( nkmer > 0 ) delete[] kmers;
    if ( nread > 0 ) delete[] reads;    
    if ( nread > 0 ) delete[] inits;
    if ( cins > 0  ) delete[] inss;
    if ( cdel > 0  ) delete[] dels;
    delete vprof;
}

void SpaPath::resetIndels()
{
    Mismatch *oins = inss;
    Mismatch *odel = dels;
    if ( cins > 0 ) delete[] oins;
    if ( cdel > 0 ) delete[] odel;
    cins = cdel = 0;
    inss = dels = NULL;
}

/**
 * k-mer existence flags.
 * 2D array - #reads by #kmers.
 * For each read, check exitence of k-mers.
 */
void SpaPath::__initFlags(bool **flags, InvertedIndex &iindex)
{
    ReadIdArray path_reads = ReadIdArray( reads, reads+nread );
    sort( path_reads.begin(), path_reads.end() );
    for ( unsigned i = 0; i < nkmer; i++ ) {
        KmerId kid = kmers[i];
        ReadIdArray crids = ReadIdArray( iindex.getValue(kid)->rid, iindex.getValue(kid)->rid + iindex.getValue(kid)->size );
        std::sort( crids.begin(), crids.end() );
        ReadIdArray comm = Set::Intersect<ReadId>( path_reads, crids, true );
        ReadIdSet scom = ReadIdSet(comm.begin(), comm.end());
		
        for ( unsigned j = 0; j < nread; j++ ) {
            if ( scom.find(reads[j]) != scom.end() ) flags[i][j] = true;
            else flags[i][j] = false;
        }
    }
}

void SpaPath::__initFlags(bool **flags, BitString *bstrs, int kmer_size)
{
    for ( unsigned i = 0; i < nkmer; i++ ) {
        KmerId kid = kmers[i];
        std::string kmer = alpha::IntegerToAminoAcid(kid, kmer_size);

        for ( unsigned j = 0; j < nread; j++ ) {
            std::string rstr = bstrs[reads[j]].toString();
            size_t found = rstr.find(kmer);
            if ( found != std::string::npos ) flags[i][j] = true;
            else flags[i][j] = false;
        }
    }
}

void SpaPath::__destroy(bool **flags)
{
    for ( unsigned i = 0; i < nkmer; i++ )
        delete[] flags[i];
    delete[] flags;

}

std::string SpaPath::__alignedString(bool **flags, unsigned rid, unsigned s, unsigned e, unsigned k)
{
    unsigned aln_len = k+(e-s);
    std::string aln_str = std::string(aln_len, '.');
    for ( unsigned j = s; j <= e; j++ ) {
        if ( flags[j][rid] ) {
            KmerType kmer = alpha::IntegerToAminoAcid(kmers[j], k);
            aln_str.replace(j-s, k, kmer);
        } 
    }
    return aln_str;
}

uList SpaPath::__mismatch(bool **flags, unsigned rid, unsigned s, unsigned e, unsigned k)
{
    std::string aln_str = __alignedString( flags, rid, s, e, k );
    uList locs;
    for ( unsigned j = 0; j < aln_str.size(); j++ ) 
        if ( aln_str[j] == '.' ) locs.push_back(j);
    return locs;
}

void SpaPath::__trim(uList &mis, bool verbose)
{
    if ( nread == 0 ) return;

    unsigned oread = nread;
    std::set<unsigned> dropset = std::set<unsigned>(mis.begin(), mis.end());
    nread -= dropset.size();

    if ( verbose ) std::cout << "# reads (old:" << oread << "\tnew:" << nread << "\tskip:" << oread-nread << ")\n";

    ReadId *orids = reads;
    int *oss = inits;

    if ( nread > 0 ) {
        reads = new ReadId[nread];
        inits = new int[nread];
    }
    
    unsigned curr = 0;
    for ( unsigned  i = 0; i < oread; i++ ) {
        if ( dropset.find(i) == dropset.end() ) {
            reads[curr] = orids[i];
            inits[curr] = oss[i];
            curr++;
        }
    }
    delete[] orids;
    delete[] oss;
}

unsigned SpaPath::__start( bool **flags, unsigned rid, unsigned s )
{
    while ( flags[s][rid] == false ) { 
        s++; 
        if ( s == nkmer ) break; // Shouldn't happen. Safety check;
    }
    return s;
}

unsigned SpaPath::__start( bool **flags, unsigned rid )
{
    unsigned s = 0;
//     while ( flags[s][rid] == false ) { s++; }
//     return s;
    return __start(flags, rid, s);
}

unsigned SpaPath::__nstart( bool **flags, unsigned rid, unsigned s )
{
    while ( s < nkmer && flags[s][rid] == true  ) { s++; }
    s++;
    while ( s < nkmer && flags[s][rid] == false ) { s++; }
    return s;
}

unsigned SpaPath::__end( bool **flags, unsigned rid, unsigned e )
{
    while ( flags[e][rid] == false ) { 
        e--; 
        if ( e == 0 ) break; // Safety check
    }
    return e;
}

unsigned SpaPath::__end( bool **flags, unsigned rid )
{
    unsigned e = nkmer-1;
    return __end(flags, rid, e);
}

int __scorePos( bool **flags,
                int i,
                unsigned s,
                unsigned e )
{
    int match = 1;
    int gapop = -1l;
    int gapex = -1;
    int score = 0;
    bool gflag = false;
    for ( unsigned a = s; a <= e; a++ ) {
//         if ( flags[a][i] == true )
//             std::cout << "1";
//         else
//             std::cout << "0";
        if ( flags[a][i] == true ) {
            gflag = false;
            score += match;
        }
        else {
            if ( gflag ) score += gapex;
            else {
                score += gapop;
                gflag = true;
            }
        }
    }
//     std::cout << "\n";
    return score;
}

IntPair SpaPath::__getRange( bool **flags,
                            unsigned i,
                            unsigned seq_len,
                            int k )
{
    IntPair se;
    
    std::vector<int> svec, evec;
    unsigned s = __start(flags, i);
    //unsigned is = s;
    unsigned e = __end(flags, i); 

    return IntPair(s,e);
    /*
    svec.push_back(s);
    evec.push_back(e);

    while( s < e ) {
        while( flags[s][i] == true && s < e ) 
            s++;
        if ( s == e ) break; // consecutive trues
   
        while ( flags[s][i] == false )
            s++;

        svec.push_back(s);
    }

    while( e > is ) {
        while( flags[e][i] == true && e > is ) 
            e--;
        if ( e == is ) break; // consecutive trues
   
        while ( flags[e][i] == false )
            e--;

        evec.push_back(e);
    }

//     std::cout << "Starts:";
//     for ( int l = 0; l < svec.size(); l++ ) 
//         std::cout << svec[l] << " ";
//     std::cout << "\n";
//     std::cout << "Ends:";
//     for ( int l = 0; l < evec.size(); l++ ) 
//         std::cout << evec[l] << " ";
//     std::cout << "\n";

    int maxscore = -1000;
    IntPair maxpair = IntPair(svec[0], evec[0]);
    for ( size_t l = 0; l < svec.size(); l++ ) {
        for ( size_t j = 0; j < evec.size(); j++ ) {
            if ( k+evec[j]-svec[l] == (int)seq_len) return IntPair(svec[l], evec[j]);
            int score = __scorePos(flags, i, svec[l], evec[j]);
            //std::cout << svec[l] << "\t" << evec[j] << "\tscore:" << score << "\n";
            if ( score > maxscore ) {
                maxscore = score;
                maxpair = IntPair(svec[l], evec[j]);
            }
        }
    }

//     std::cout << "Max score:" << maxscore << "\n";
//     std::cout << "Max pair: " << maxpair.first << " " << maxpair.second << "\n";

    return maxpair;
//     if ( k+(e-s) == seq_len ) return IntPair(s,e);

//     //if ( s == e ) return IntPair(s,e);
//     int max = __scorePos(flags, i, s, e);
    
//     while(1) {
//         if ( s >= e ) return IntPair(e, e);
//         if ( flags[s+1][i] == true ) break;
//         s = __start(flags, i, s+1);
//     }
//     while(1) {
//         if ( e <= s ) return IntPair(s, s);
//         if ( flags[e-1][i] == true ) break;
//         e = __end(flags, i, e-1);
//     }
//     return IntPair(s,e);
*/
}

void SpaPath::validate(BitString *bstrs, Param &param)
{
    if ( nread == 0 ) return;

    int mstart = inits[0];
    for ( size_t i = 1; i < nread; i++ ) 
        if ( inits[i] < mstart ) mstart = inits[i];

    if ( mstart > 0 ) {
        if ( param.verbose ) std::cout << "MStart:" << mstart << "\n";
        
        KmerId *okmers = kmers;
        nkmer -= mstart;
        kmers = new KmerId[nkmer];
        for ( size_t j = 0; j < nkmer; j++ )
            kmers[j] = okmers[mstart+j];
        delete[] okmers;

        for ( size_t i = 0; i < nread; i++ ) {
            inits[i] -= mstart;
        }
        for ( size_t i = 0; i < cins; i++ ) {
            inss[i].rpos -= mstart;
        }
        for ( size_t i = 0; i < cdel; i++ ) {
            dels[i].rpos -= mstart;
        }
    }
    
    int rlen = nkmer + param.kmer_size - 1;
    int mend = inits[0] + bstrs[reads[0]].getSize();
    for ( size_t i = 1; i < nread; i++ ) {
        int clen = inits[i] + bstrs[reads[i]].getSize();
        if ( mend < clen ) mend = clen;
    }

    if ( rlen > mend ) {
        if ( param.verbose ) std::cout << "MEnd:" << rlen - (rlen-mend) << "\n";
        int diff = rlen - mend;
        KmerId *okmers = kmers;
        nkmer -= diff;        
        kmers = new KmerId[nkmer];
        for ( size_t j = 0; j < nkmer; j++ )
            kmers[j] = okmers[j];
        delete[] okmers;
    }
}

//==============================================================================
// To do. 
// Determine type of mismatch here.
//==============================================================================
//void SpaPath::prealign(InvertedIndex &iindex, unsigned k, char **seqs, double pcut, bool verbose)
void SpaPath::align(InvertedIndex &iindex, BitString *bstrs, Param &param)
{
    /* 2D kmer existence flags */
    bool **flags;
    flags = new bool *[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ ) flags[i] = new bool[nread];
    __initFlags(flags, iindex);
    //__initFlags(flags, bstrs, param.kmer_size);

    std::string ref = biostr::getSequenceString(kmers, nkmer, param.kmer_size);
    uList mis;
    std::list<Mismatch>ilist, dlist; // temprary lists;
    if ( param.verbose ) std::cout << "# reads:" << nread << "\n";
    for ( unsigned i = 0; i < nread; i++ ) {
        AlignSummary summary;
        bool success = __alignReadToPath(summary, ref, i, flags, bstrs, ilist, dlist, param);
        if ( !success ) mis.push_back(i); 
    }

    resetIndels();
    __listToArray(ilist, INSERTION);
    __listToArray(dlist, DELETION);
    __trim(mis, param.verbose);
    __destroy(flags);
}  

void SpaPath::align(BitString *bstrs, Param &param)
{
    /* 2D kmer existence flags */
    bool **flags;
    flags = new bool *[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ ) flags[i] = new bool[nread];
    __initFlags(flags, bstrs, param.kmer_size);

    std::string ref = biostr::getSequenceString(kmers, nkmer, param.kmer_size);
    uList mis;
    std::list<Mismatch>ilist, dlist; // temprary lists;
    if ( param.verbose ) std::cout << "# reads:" << nread << "\n";
    for ( unsigned i = 0; i < nread; i++ ) {
        AlignSummary summary;
        bool success = __alignReadToPath(summary, ref, i, flags, bstrs, ilist, dlist, param);
        if ( !success ) mis.push_back(i); 
    }

    resetIndels();
    __listToArray(ilist, INSERTION);
    __listToArray(dlist, DELETION);
    __trim(mis, param.verbose);
    __destroy(flags);
}

bool SpaPath::__short( unsigned seq_len, unsigned aln_len, Param &param )
{
    if ( aln_len < unsigned(param.read_align_ratio * seq_len) ) {
        if ( param.verbose ) std::cout << "\tshort\n"; 
        return true;
    }
    return false;
}

bool SpaPath::__weak(IntPair &se, unsigned seq_len, unsigned aln_len, bool **flags, unsigned i, Param &param)
{
    uList mmp =__mismatch(flags, i, se.first, se.second, param.kmer_size);
    if ( (aln_len-mmp.size()) < unsigned(param.read_align_score*aln_len) ) {
        if ( param.verbose ) std::cout << "\tmultiple mismatches:" << mmp.size() << "\n";
        return true;
    }
    return false;
}

AlignSummary SpaPath::__compareBase(std::string &query, std::string &sbjct)
{
    AlignSummary summary;
    summary.length = query.length();
    summary.positive = scoring::countPositive(query, sbjct, BLOSUM62);
    summary.posrate = (double)summary.positive/summary.length;
    return summary;
}

bool SpaPath::__alignReadToPath(AlignSummary &summary, std::string &ref, unsigned i, bool **flags, BitString *bstrs, std::list<Mismatch> &ilist, std::list<Mismatch> &dlist, Param &param)
{
    std::string query = bstrs[reads[i]].toString();
    unsigned     qlen = query.size();
    IntPair        se = __getRange(flags, i, qlen, param.kmer_size);
    unsigned        s = se.first;
    unsigned        e = se.second;
    unsigned     aln_len = param.kmer_size+(e-s);

    if ( param.verbose ) std::cout << i << "\tr:" << reads[i] << "\tl:" << qlen << "\ts:" << s << "\te:" << e << "\ta:" << aln_len << "\n";

    if ( e < s ) return false; // Safety check
    if ( __short( qlen, aln_len, param ) ) return false;
    if ( aln_len < qlen && __weak( se, qlen, aln_len, flags, i, param ) ) return false;

    /**
     * Initialize alignment start position of the current read
     */
    inits[i] = s;


    /**
     * First, check substring match 
     */
    size_t pos = ref.find(query);
    if ( pos != std::string::npos ) {
        if ( param.verbose ) std::cout << "\tSubstring match to sbjct at:" << pos << "\n";
        inits[i] = pos;
        return true;
    }

    /** 
     * Check substring match
     * query length > aligned region length
     */
    if ( aln_len < qlen ) {
        assert( s+aln_len <= ref.size() );
        assert(s >= 0);
        std::string sbjct = ref.substr(s, aln_len);
        size_t pos = query.find(sbjct);
        if ( pos != std::string::npos ) {
            if ( param.verbose ) std::cout << "\tSubstring match to query at:" << inits[i] - (int)pos << "\n";
            inits[i] = inits[i] - (int)pos;
            return true;
        }
    }
    
    

    /**
     * Compare by substitution 
     */
    if ( aln_len == qlen ) {
        assert( s+aln_len <= ref.size() );
        assert(s >= 0);
        std::string sbjct = ref.substr(s, aln_len);
        summary = __compareBase(query, sbjct);
        if ( summary.posrate < param.read_align_score ) return false;
        return true;
    }
    
    /** 
     * Try to compare base first  
     */
    while(true) {
        aln_len = param.kmer_size+(e-s);
        if ( aln_len < qlen ) break;
        
        //-----------------------------
        // Same length of query & sbjct
        //-----------------------------
        assert( s+qlen <= ref.size() );
        assert(s >= 0);
        std::string sbjct = ref.substr(s, qlen);
        summary = __compareBase(query, sbjct);
        if ( summary.posrate >= param.read_align_score ) {
            if ( param.verbose ) std::cout << "\tSubstitution success at:" << s << "\tscore:" << summary.posrate << "\n";
            inits[i] = s; 
            return true;
        }
        s = __nstart(flags, i, s );
        if ( s == nkmer ) break;
    }
    
    
    /** 
     * Now align anyway 
     */
    int init = se.first;
    aln_len = param.kmer_size+(se.second-se.first);
    assert( init+aln_len <= ref.size() );
    assert( init >= 0 );
    std::string sbjct = ref.substr(init, aln_len);
    bool good = __doAlignment(summary, query, sbjct, param);
    if ( !good ) return false;
    //inits[i] = 
    __adjustReadStart( i, init, summary, param.verbose );
    __updateIndel(i, init, summary, ilist, dlist);
            
    return true;

//     if ( aln_len < seq_len ) {
//         if ( aln_len < unsigned(pcut * seq_len) ) {
//             if ( verbose ) std::cout << "\tshort\n"; 
//             mis.push_back(i); continue;
//         } else {
//             uList mmp =__mismatch(flags, i, s, e, k);
//             if ( (aln_len-mmp.size()) < unsigned(pcut*aln_len) ) {
//                 if ( verbose ) std::cout << "\tmultiple mismatches:" << mmp.size() << "\n";
//                 if ( verbose ) {
//                     //std::cout << "\tseq:" << seqs[reads[i]] << "\n";
//                     std::cout << "\tseq:" << seq << "\n";
//                     std::cout << "\tref:" << ref.substr(s, aln_len) << "\n";
//                 }
//                 mis.push_back(i); continue;
//             }
//         }
//     }

//     // no alignment
//     else if ( aln_len == seq_len )  { 
//         std::string rstr = ref.substr(s, aln_len);
//         //string qstr = string(seqs[reads[i]]);
//         std::string qstr = seq;
//         if ( rstr != qstr ) {
//             int nmatch = 0;
//             for ( size_t p = 0; p < aln_len; p++ )
//                 if ( rstr[p] == qstr[p] ) nmatch++;
//             if ( nmatch < aln_len*pcut ) {
//                 mis.push_back(i);
//                 if ( verbose ) std::cout << "\tmultiple mismatches\n";
//                 if ( verbose ) {
//                     //cout << "\tseq:" << seqs[reads[i]] << "\n";
//                     std::cout << "\tseq:" << qstr << "\n";
//                     std::cout << "\tref:" << ref.substr(s, aln_len) << "\n";
//                 }
//             }
//         }
//         continue;
//     }

//     /*
//       L-----------EQQQANAEPEQQVDPRKAAVEAAIARAKARKL------
//       ||| |||||| |||||||||||||||||||||
//       -AAIARAKARKREQQPANAEPEEQVDPRKAAVEAAIARAKARKLEQQQAN
//     */

//     //-------
//     // indels
//     //-------
//     //string query = std::string(seqs[reads[i]]);
//     std::string query = seq;
//     std::string match = ref.substr(inits[i], aln_len);
//     //GlobalAlignPair paln = GlobalAlignPair(query, match);
//     GlobalAlignPair paln = GlobalAlignPair(match, query);
//     AlignSummary summary = paln.getSummary();

//     size_t aln_beg = summary.range.first;
//     size_t aln_end = summary.range.second;
//     aln_len = aln_end - aln_beg + 1;
                    
//     if ( verbose ) {
//         std::cout << paln.getAlignment();
            
//         std::cout << "\t#ins:" << summary.ins.size() << "\t";
//         std::cout << "\t#del:" << summary.del.size() << "\t";
//         std::cout << "\t#mat:" << summary.match << "\t";
//         std::cout << "\t#pos:" << summary.positive << "\t";
//         std::cout << "\t#mms:" << summary.mismatch << "\n";
//         std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//         std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
//     }


//     if ( aln_len < unsigned(pcut * query.size()) ) {
//         if ( verbose ) std::cout << "\talignment - short\n";
//         mis.push_back(i);
//         continue;
//     }
//     else if ( aln_len * pcut  >  query.size() ) {
//         if ( verbose ) std::cout << "\talignment - long\n";
//         mis.push_back(i);
//         continue;
//     }
        
//     //if ( summary.positive < aln_len * pcut ) {
//     if ( (double)summary.positive <  aln_len * pcut ) {
//         mis.push_back(i);
//         if ( verbose ) std::cout << "\talignment - weak\n";
//         continue;
//     }

//     int rpos = inits[i];
        
//     // adjust with indels
//     inits[i] = adjustReadStart( inits[i], summary, verbose );
        
//     if ( verbose ) {
//         if ( summary.ilist.size() > 0 ) {
//             std::cout << "\tInsertions:\t";
//             for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
//                 //std::cout << "aln:" << it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << it->ref_pos + inits[i] << "\t";
//                 std::cout << "aln:" << it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << rpos + it->ref_pos << "\t";
//             }
//         } 
//         if ( summary.dlist.size() > 0 ) {
//             std::cout << "\tDeletions:\t";
//             for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
//                 //std::cout << "aln:" <<it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << it->ref_pos + inits[i] << "\t";
//                 std::cout << "aln:" <<it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << rpos + it->ref_pos << "\t";
//             }
//         }
//         std::cout << "\n";
//     }


}

void SpaPath::__updateIndel(unsigned i,
                           int init,
                           AlignSummary &summary,
                           std::list<Mismatch> &ilist,
                           std::list<Mismatch> &dlist)
{
    AlignPosList::iterator it;
    for ( it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
        ilist.push_back( Mismatch(reads[i], it->seq_pos, init + it->ref_pos) );
    }
    for ( it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
        dlist.push_back( Mismatch(reads[i], it->seq_pos, init + it->ref_pos) );
    }
}

//int 
void SpaPath::__adjustReadStart( int index,
                                int init, 
                                AlignSummary &summary,
                                bool verbose )
{
    if ( verbose ) {
        
    }
    if ( summary.lgap.first == 0 && summary.lgap.second > 0 ) {
        if ( verbose ) std::cout << "Read:" << reads[index];
        if ( verbose ) std::cout << "\tstart:" << init;
        if ( verbose ) std::cout << " -> ";
        int diff = (summary.lgap.second-summary.lgap.first);
        init = init + diff;
        if ( summary.range.first != diff ) init += summary.range.first;
        if ( verbose ) std::cout << init;
        if ( verbose ) std::cout << "\n";
    } else if ( summary.lgap.first > 0 && summary.lgap.second == 0 ) {
        if ( verbose ) std::cout << "Read:" << reads[index];
        if ( verbose ) std::cout << "\tstart:" << init;
        if ( verbose ) std::cout << " -> ";
        int diff =(summary.lgap.first-summary.lgap.second);
        init = init - diff;
        if ( summary.range.first != diff ) init += summary.range.first;
        if ( verbose ) std::cout << init;
        if ( verbose ) std::cout << "\n";
    }
    //return init;
    inits[index] = init;
}

void SpaPath::__listToArray(std::list<Mismatch> &clist, int type)
{
    int osize;
    type == INSERTION ? osize = cins : osize = cdel;

    unsigned nsize = clist.size();
    //if ( nsize == 0 ) return;

    if ( type == INSERTION ) {
        cins = nsize;
        Mismatch *oinss = inss;

        if ( cins > 0 ) {
            inss = new Mismatch[cins];
            int i = 0;
            for ( std::list<Mismatch>::iterator it = clist.begin(); it != clist.end(); ++it ) {
                inss[i] = *it; i++;
            }
        }
        else inss = NULL;
        
        if ( osize > 0 ) delete[] oinss;

//         if ( osize > 0 ) {
//             delete[] inss; inss = NULL;
//         }
//         if ( nsize == 0 ) return;

//         inss = new Mismatch[nsize];
//         int i = 0;
//         for ( std::list<Mismatch>::iterator it = clist.begin(); it != clist.end(); ++it ) {
//             inss[i] = *it; i++;
//         }
    }
    else {
        cdel = nsize;
        Mismatch *odels = dels;

        if ( cdel > 0 ) {
            dels = new Mismatch[cdel];
            int i = 0;
            for ( std::list<Mismatch>::iterator it = clist.begin(); it != clist.end(); ++it ) {
                dels[i] = *it; i++;
            }
        }
        else dels = NULL;

        if ( osize > 0 ) delete[] odels;

//         if ( osize > 0 ) {
//             delete[] dels; dels = NULL;
//         }
//         if ( nsize == 0 ) return;

//         dels = new Mismatch[nsize];
//         int i = 0;
//         for ( std::list<Mismatch>::iterator it = clist.begin(); it != clist.end(); ++it ) {
//             dels[i] = *it; i++;
//         }
    }
}

// void SpaPath::__updateindel(uList &mmp, unsigned i )
// {
    
// }

// //void SpaPath::__alnStartPos( char** seqs, unsigned i, unsigned s, unsigned k )
// //unsigned SpaPath::__alnStartPos( char** seqs, unsigned i, unsigned s, unsigned k, bool verbose )
// unsigned SpaPath::__alnStartPos( BitString *bstrs, unsigned i, unsigned s, unsigned k, bool verbose )
// {
//     //inits[i] = s;
//     KmerType skmer = alpha::IntegerToAminoAcid(kmers[s], k);
//     //string seq = std::string(seqs[reads[i]]);
//     std::string seq = bstrs[reads[i]].toString();
    
//     if ( verbose ) std::cout << "s:" << s << "\n";
//     //unsigned os = inits[i];
//     unsigned os = s;
//     unsigned bf =  __fitstart(skmer, seq, verbose);
//     if ( verbose ) std::cout << "bf:" << bf << "\n";
//     //inits[i] -= bf;
//     if ( verbose ) std::cout << i << "\tos:" << os << "\tbf:" << bf << "\tns:" << os-bf << "\n";
//     if ( bf > 0 ) {
//         if ( verbose ) std::cout << (int)s-(int)bf << "\n";
//         if ( (int)s-(int)bf < 0 ) {
//             bf = s; 
//             if ( verbose ) std::cout << "new bf:" << bf << "\n"; 
//         }
//         if ( verbose ) std::cout << "\tolds:" << skmer << "\n";
//         KmerType nkmer = alpha::IntegerToAminoAcid(kmers[s-bf], k);
//         if ( verbose ) std::cout << "\tnews:" << nkmer << "\n";
//         if ( verbose ) std::cout << "\tread:" << seq << "\n";
//     }
//     return bf;
// }

bool SpaPath::__doAlignment(AlignSummary &summary, std::string &query, std::string &sbjct, Param &param)
{
//     std::string query = seq;
//     std::string match = ref.substr(inits[i], aln_len);
//     //GlobalAlignPair paln = GlobalAlignPair(query, match);
    GlobalAlignPair paln = GlobalAlignPair(sbjct, query);
    //AlignSummary 
    summary = paln.getSummary();

    size_t aln_beg = summary.range.first;
    size_t aln_end = summary.range.second;
    size_t aln_len = aln_end - aln_beg + 1;
                    
    if ( param.verbose ) {
        std::cout << paln.getAlignment();
            
        std::cout << "\t#ins:" << summary.ins.size() << "\t";
        std::cout << "\t#del:" << summary.del.size() << "\t";
        std::cout << "\t#mat:" << summary.match << "\t";
        std::cout << "\t#pos:" << summary.positive << "\t";
        std::cout << "\t#mms:" << summary.mismatch << "\n";
        std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
        std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
    }


    if ( aln_len < unsigned(param.read_align_ratio * query.size()) ) {
        if ( param.verbose ) std::cout << "\talignment - short\n";
        return false;
    }
    else if ( aln_len * param.read_align_ratio  >  query.size() ) {
        if ( param.verbose ) std::cout << "\talignment - long\n";
        return false;
    }

    if ( (double)summary.positive <  aln_len * param.read_align_score ) {
        if ( param.verbose ) std::cout << "\talignment - weak\n";
        return false;
    }

    return true;

}

int SpaPath::__matchlength( std::string &str1, std::string &str2 )
{
    int len = 0;
    for ( size_t i = 0; i < str1.size(); i++ )
        if ( str1[i] == str2[i] ) len++;
    return len;
}


void SpaPath::dump(std::ostream &out)
{
    unsigned len = (unsigned)strlen(consensus);
    out.write((char*)&len, sizeof(unsigned));
    char ch;
    for ( unsigned i = 0; i < len; i++ ) {
        ch = consensus[i];
        out.write((char*)&ch, sizeof(char));
    }
    
    out.write((char*)&nkmer, sizeof(unsigned));
    out.write((char*)&nread, sizeof(unsigned));
    out.write((char*)&cins,  sizeof(unsigned));
    out.write((char*)&cdel,  sizeof(unsigned));
    
    for ( unsigned i = 0; i < nkmer; i++ )
        out.write((char*)&kmers[i], sizeof(KmerId));
    for ( unsigned i = 0; i < nread; i++ )
        out.write((char*)&reads[i], sizeof(ReadId));
    for ( unsigned i = 0; i < nread; i++ )
        out.write((char*)&inits[i], sizeof(int));

    for ( unsigned i = 0; i < cins; i++ )
        inss[i].dump(out);
    for ( unsigned i = 0; i < cdel; i++ )
        dels[i].dump(out);

    vprof->dump(out);
}

void SpaPath::load(std::istream &in)
{
    unsigned len;
    in.read((char*)&len, sizeof(unsigned));
    //std::cout << "length:" << len << "\n";
    consensus = new char[len+1];
    //in.read((char*)&consensus, sizeof(char*));
    //in.read((char*)&consensus, sizeof(char*)*(len+1));
    //std::cout << "Consensus:" << consensus << "\n";
    char ch;
    for ( unsigned i = 0; i < len; i++ ) {
        //in.read((char*)&consensus[i], sizeof(char));
        in.read((char*)&ch, sizeof(char));
        consensus[i] = ch;
    }
    consensus[len] = '\0';

    //std::cout << "Consensus:" << consensus << "\n";
    in.read((char*)&nkmer, sizeof(unsigned));
    in.read((char*)&nread, sizeof(unsigned));
    in.read((char*)&cins,  sizeof(unsigned));
    in.read((char*)&cdel,  sizeof(unsigned));

    
    if ( nkmer > 0 ) kmers = new KmerId[nkmer];
    if ( nread > 0 ) reads = new ReadId[nread];
    if ( nread > 0 ) inits = new int[nread];

    for ( unsigned i = 0; i < nkmer; i++ )
        in.read((char*)&kmers[i], sizeof(KmerId));
    for ( unsigned i = 0; i < nread; i++ )
        in.read((char*)&reads[i], sizeof(ReadId));
    for ( unsigned i = 0; i < nread; i++ )
        in.read((char*)&inits[i], sizeof(int));


    if ( cins > 0 ) inss = new Mismatch[cins];
    for ( unsigned i = 0; i < cins; i++ )
        inss[i].load(in);

    if ( cdel > 0 ) dels = new Mismatch[cdel];
    for ( unsigned i = 0; i < cdel; i++ )
        dels[i].load(in);

    vprof = new ProfileVector();
    vprof->load(in);
    
}

void SpaPath::prepend(SpaPath *oaln, int lgap, int kmer_size)
{
    KmerId *okmers = kmers;
    std::string this_seq = biostr::getSequenceString(kmers, nkmer, kmer_size);
    std::string that_seq = biostr::getSequenceString( oaln->getKmers(), oaln->getKmerCount(), kmer_size );
    
    this_seq = that_seq + this_seq;
    KmerArray kids = biostr::getKmers(this_seq, kmer_size);
    nkmer = (unsigned)kids.size();
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ )
        kmers[i] = kids[i];
    delete[] okmers;

    updateConsensusSequence( this_seq.c_str() );


//     KmerId *okmers = kmers;
//     nkmer += lgap;
//     kmers = new KmerId[nkmer];
    
//     KmerId *others = oaln->getKmers();
//     for ( int j = 0; j < lgap; j++ ) kmers[j] = others[j];
//     for ( int j = lgap; j < (int)nkmer; j++ ) kmers[j] = okmers[j-lgap];
    
//     delete[] okmers;

//     std::string nstr = biostr::getSequenceString( kmers, nkmer, kmer_size );
//     updateConsensusSequence( nstr.c_str() );
}

void SpaPath::prepend(std::string &that_seq, size_t kmer_size)
{
    KmerId *okmers = kmers;
    std::string this_seq = biostr::getSequenceString(kmers, nkmer, kmer_size);
    this_seq = that_seq + this_seq;
    KmerArray kids = biostr::getKmers(this_seq, kmer_size);
    nkmer = (unsigned)kids.size();
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ )
        kmers[i] = kids[i];
    delete[] okmers;

    updateConsensusSequence( this_seq.c_str() );
}


// KmerArray SpaPath::getKmerArray( std::string &seq, int k )
// {
//     std::list<KmerId> kids;
//     for (size_t i = 0; i <= seq.size() - k; i++) {
//         KmerType kmer = seq.substr(i,k);
//         kids.push_back( alpha::AminoAcidToInteger<KmerId>(kmer) );
//     }
//     return KmerArray( kids.begin(), kids.end() );
// }

void SpaPath::append(SpaPath *oaln, int tgap, int kmer_size)
{
     KmerId *okmers = kmers;
//     size_t onkmer = nkmer;
    
//     nkmer += tgap;
//     kmers = new KmerId[nkmer];
    
    std::string this_seq = biostr::getSequenceString(kmers, nkmer, kmer_size);
    std::string that_seq = biostr::getSequenceString(oaln->getKmers(), oaln->getKmerCount(), kmer_size);

    assert( tgap <= (int)that_seq.size() );
    this_seq += that_seq.substr( that_seq.length()-tgap, tgap );
    KmerArray kids = biostr::getKmers(this_seq, kmer_size);
    nkmer = (unsigned)kids.size();
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ )
        kmers[i] = kids[i];

//     nkmer = this_seq.length()-kmer_size+1;
    
//     for ( size_t j = 0; j < onkmer; j++ ) kmers[j] = okmers[j];
//     KmerId *others = oaln->getKmers();
//     size_t other_nk = oaln->getKmerCount();
//     for ( int j = 0; j < tgap; j++ ) {
//         kmers[onkmer+j] = others[other_nk-tgap+j];
//         std::cout << "\t" << j << ":adding kmer:" << 
//             alpha::IntegerToAminoAcid<KmerId>(others[other_nk-tgap+j], 6) << "\n";
//    }
    delete[] okmers;

    //std::string nstr = biostr::getSequenceString( kmers, nkmer, kmer_size );
    updateConsensusSequence( this_seq.c_str() );//nstr.c_str() );
}

void SpaPath::append(std::string &that_seq, size_t kmer_size)
{
    KmerId *okmers = kmers;
    std::string this_seq = biostr::getSequenceString(kmers, nkmer, kmer_size);
    this_seq += that_seq;
    KmerArray kids = biostr::getKmers(this_seq, kmer_size);
    nkmer = (unsigned)kids.size();
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ )
        kmers[i] = kids[i];
    delete[] okmers;

    updateConsensusSequence( this_seq.c_str() );
}

// void SpaPath::append(KmerId *okmers, size_t onkmer)
// {
//     KmerId *old_kmers = kmers;
//     size_t old_nkmer  = nkmer;

//     nkmer += onkmer;
//     kmers = new KmerId[nkmer];
    
//     for ( size_t j = 0; j < old_nkmer; j++ ) kmers[j] = old_kmers[j];
//     for ( size_t j = 0; j < onkmer; j++ ) kmers[old_nkmer+j] = okmers[j];
    
//     delete[] okmers;
// }

void SpaPath::mergeReads(SpaPath *oaln)
{
    size_t old_nreads = nread;
    ReadId *old_reads = reads;
    int *old_inits = inits;

    size_t other_nreads = oaln->getReadCount();
    ReadId *other_reads = oaln->getReads();
    int *other_inits = oaln->getInits();

    nread += other_nreads;
    reads = new ReadId[nread];
    inits = new int[nread];
    for ( size_t i = 0; i < old_nreads; i++ ){
        reads[i] = old_reads[i];
        inits[i] = old_inits[i];
    }
    for ( size_t i = 0; i < other_nreads; i++ ){
        reads[i+old_nreads] = other_reads[i];
        inits[i+old_nreads] = other_inits[i];
        //std::cout << "New member:" << other_reads[i] << "\t" << other_inits[i] << "\n";
    }
    
    delete[] old_reads;
    delete[] old_inits;
}

void SpaPath::mergeIndels(SpaPath *oaln)
{
    std::list<Mismatch> ilist, dlist;
//     for (int i = 0; i < cins; i++) ilist.push_back(inss[i]);
//     for (int i = 0; i < cdel; i++) dlist.push_back(dels[i]);
    if (cins) ilist = std::list<Mismatch>(inss, inss+cins);
    if (cdel) dlist = std::list<Mismatch>(dels, dels+cdel);
    
    Mismatch *oins = oaln->getInsertions();
    Mismatch *odel = oaln->getDeletions();
    unsigned ocins = oaln->countInsertions();
    unsigned ocdel = oaln->countDeletions();
    for ( unsigned i = 0; i < ocins; i++ )
        ilist.push_back(oins[i]);
    for ( unsigned i = 0; i < ocdel; i++ )
        dlist.push_back(odel[i]);
    
    // Mismatch *oins = inss;
    //Mismatch *odel = dels;
    //oins = odel = NULL;

    resetIndels();
//     if ( cins > 0 ) delete[] inss;
//     if ( cdel > 0 ) delete[] dels;
//     cins = cdel = 0;
//     inss = dels = NULL;

    __listToArray(ilist, INSERTION);
    __listToArray(dlist, DELETION);

}

void SpaPath::adjustInit(int i, int g)
{
    inits[i] -= g;
}

void SpaPath::adjustStartPositions(int offset)
{
    for ( unsigned i = 0; i < nread; i++ ) {
        inits[i] += offset;
    }
}

void SpaPath::adjustPositions(int offset)
{
    //std::cout << "adjusting to " << offset << "\n";
    for ( unsigned i = 0; i < nread; i++ ) {
        //std::cout << "Read:" << reads[i] << "\tOld:" << inits[i] << "\t";
        inits[i] += offset;
        //std::cout << "new:" << inits[i] << "\n";
    }
    for ( unsigned i = 0; i < cins; i++ )
        inss[i].rpos += offset;
    for ( unsigned i = 0; i < cdel; i++ )
        dels[i].rpos += offset;
}

//
// delete insertion in reference if the position is bigger than ref size.
//
void SpaPath::trimIndels(size_t spos)
{
    std::list<Mismatch>ilist, dlist; // temprary lists;

    for ( unsigned i = 0; i < cins; i++ ) {
        if ( inss[i].rpos < (int)spos ) ilist.push_back(inss[i]);
    }
    for ( unsigned i = 0; i < cdel; i++ ) {
        if ( dels[i].rpos < (int)spos ) dlist.push_back(dels[i]);
    }

    Mismatch *oins = inss;
    Mismatch *odel = dels;
    oins = odel = NULL;
    if ( cins > 0 ) delete[] oins;
    if ( cdel > 0 ) delete[] odel;
    cins = cdel = 0;

    __listToArray(ilist, INSERTION);
    __listToArray(dlist, DELETION);
}

AlignSummary SpaPath::__realignRead( std::string &query, int start , int end, bool verbose)
{
    assert(start >= 0);
    assert(end>=start);
    assert(strlen(consensus)>0);

    //int last = query.size()*2 - start;
    //if ( last >= strlen(consensus) ) last = strlen(consensus)-1;
    std::string sbjct = std::string(consensus).substr(start, end-start+1);


    GlobalAlignPair paln = GlobalAlignPair(sbjct, query);
    AlignSummary summary = paln.getSummary();
    //size_t aln_beg = summary.range.first;
    //size_t aln_end = summary.range.second;
    //size_t aln_len = aln_end - aln_beg + 1;
                    
    //if ( param.verbose ) {
    if ( verbose ) {
        std::cout << paln.getAlignment();
            
        std::cout << "\t#ins:" << summary.ins.size() << "\t";
        std::cout << "\t#del:" << summary.del.size() << "\t";
        std::cout << "\t#mat:" << summary.match << "\t";
        std::cout << "\t#pos:" << summary.positive << "\t";
        std::cout << "\t#mms:" << summary.mismatch << "\n";
        std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
        std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
        //}
    }
        //__updateIndel(i, start, summary, ilist, dlist);

        return summary;
}

// void SpaPath::verifyGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type )
// {
//     std::cout << "Verifying carry-over gaps\n";

//     std::tr1::unordered_map<ReadId, std::list<int> > delmap, insmap;
//     makeMismatchMap(delmap, DELETION);
//     makeMismatchMap(insmap, INSERTION);

//     std::list<Mismatch> dlist, ilist;

//     for ( size_t i = range.first; i < range.second; i++ ) {
//         std::string query = bstrs[reads[i]].toString();
//         unsigned     qlen = query.size();
        
//         if ( insmap.find(reads[i]) != insmap.end() ||
//              delmap.find(reads[i]) != delmap.end() ) {
//             if ( insmap.find(reads[i]) != insmap.end() ) 
//                 std::cout << "Indel conflict:" << reads[i] << "\n";
//             if ( delmap.find(reads[i]) != delmap.end() ) 
//                 std::cout << "Deldel check:" << reads[i] << "\n";
            
//             dropMismatchRead(i);
// //             int s = inits[i]-qlen; 
// //             if ( s < 0 ) s = 0;
// //             int e = inits[i]+2*qlen;
// //             if ( e >= strlen(consensus) ) e =  strlen(consensus)-1;
//             int s = 0;
//             int e = strlen(consensus)-1;
//             AlignSummary summary = __realignRead(query, s, e);
//             //__updateIndel(i, inits[i], summary, ilist, dlist);
//             __updateIndel(i, inits[i], summary, ilist, dlist);
//         }
//     }

//     /* Is it necessary to put in order??? */
//     /* If so, sort them */
//     for ( int j = 0; j < cdel; j++ ) 
//         dlist.push_front(dels[j]);

//     for ( int j = 0; j < cins; j++ ) 
//         ilist.push_front(inss[j]);
    
//     if ( cdel > 0 ) {
//         delete[] dels; 
//         cdel = 0;
//     }

//     if ( cins > 0 ) {
//         delete[] inss; 
//         cins = 0;
//     }

//     __listToArray(dlist, DELETION);
//     __listToArray(ilist, INSERTION);

// }

void SpaPath::updateInsertionPos( AlignPosList &nlist )
{
    nlist.clear();

    std::string cstr = std::string(consensus);
    for ( size_t i = 0; i < cstr.size(); i++ ) {
        if ( cstr[i] == '-' ) 
        //if ( cstr[i] == 'X' ) 
            nlist.push_back( AlignIndex(i, i, i) );
    }
}

void SpaPath::makeGappyConsensus(AlignPosList &ninss, Param &param)
{
    std::string nstr = biostr::getSequenceString( kmers, nkmer, param.kmer_size );
    if ( param.verbose ) std::cout << "Consensus:" << nstr << "\n";
    for ( AlignPosList::reverse_iterator it = ninss.rbegin(); it != ninss.rend(); ++it ) {
        assert(it->ref_pos < (int)nstr.size());
        nstr.insert( it->ref_pos, "-" );  
        //nstr.insert( it->ref_pos, "X" );  
        //if ( param.verbose ) std::cout << it->ref_pos << "\n";  
    }
    if ( param.verbose) std::cout << "New Consensus:" << nstr << "\n";
    KmerArray nkmers = biostr::getKmers(nstr, param.kmer_size);
    updateKmers( &nkmers[0], nkmers.size() );
    updateConsensusSequence( nstr.c_str() );
}

bool SpaPath::__badAlignment( AlignSummary &summary, int qlen, Param &param )
{
    std::map<int, bool> gap_map;
    for ( size_t i = 0; i < strlen(consensus); i++ ) {
        if ( consensus[i] == '-' ) gap_map.insert( std::pair<int, bool>(i, true) );
        else gap_map.insert( std::pair<int, bool>(i, false) );
    }

    for ( AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) 
        if ( gap_map[ it->ref_pos ] ) {
            if ( param.verbose ) std::cout << "Gaps in another gap\n";
            return true;
        }

    int ngap = 0;
    for ( int i = summary.range.first; i <= summary.range.second; i++ )
        if ( gap_map[i] ) ngap++;

    int nlen = qlen - ngap;
    //assert(nlen>0);
    if (nlen<=0) return false;
    double score = (double)summary.positive/nlen;
    if ( score  <  param.read_align_score ) {
        if ( param.verbose ) std::cout << "\talignment - weak:\tqlen:" << qlen << "\tnlen:" << nlen << "\tscore:" << score << "\n";
        return true;
    }

    return false;
}


std::list<int> SpaPath::updateGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type, Param &param )
{
    std::list<int> bad_index;

    // increment start
    for ( int i = range.first; i < range.second; i++ ) 
        for ( AlignPosList::iterator it = nlist.begin(); it != nlist.end(); ++it ) 
            if ( inits[i] >= it->ref_pos ) 
                inits[i]++; 
    
    if ( type == INSERTION ) {
        updateInsertionPos( nlist );
    }

    std::tr1::unordered_map<ReadId, std::list<int> > delmap, insmap;
    makeMismatchMap(delmap, DELETION);
    makeMismatchMap(insmap, INSERTION);

    std::list<Mismatch> dlist, ilist;
    for ( int i = range.first; i < range.second; i++ ) {
        std::string query = bstrs[reads[i]].toString();
        unsigned     qlen = query.size();

        if ( insmap.find(reads[i]) != insmap.end() || delmap.find(reads[i]) != delmap.end() ) {
            if ( insmap.find(reads[i]) != insmap.end() ) {
                if (param.verbose) {
                    std::cout << "Indel conflict:" << reads[i] << "\n";
                    for ( std::list<int>::iterator lt = insmap[reads[i]].begin(); lt != insmap[reads[i]].end(); ++lt )
                        std::cout << *lt << "\t";
                    std::cout << "\n";
                }
            }
            if ( delmap.find(reads[i]) != delmap.end() ) {
                if ( param.verbose ) {
                    std::cout << "Deldel check:" << reads[i] << "\n";
                    for ( std::list<int>::iterator lt = delmap[reads[i]].begin(); lt != delmap[reads[i]].end(); ++lt )
                        std::cout << *lt << "\t";
                    std::cout << "\n";
                }
            }
            
            dropMismatchRead(i);
            int s = 0;
            int e = strlen(consensus)-1;
            AlignSummary summary = __realignRead(query, s, e, param.verbose);
            if ( __badAlignment(summary, qlen, param) ) {
                //AlignSummary summary;
                //std::string sbjct = std::string(consensus);
                //bool good = __doAlignment(summary, query, sbjct, param);
                //if ( !good ) {
                bad_index.push_back(i);
                continue;
            }
            inits[i] = summary.range.first;
            __updateIndel(i, 0, summary, ilist, dlist);
            continue;
        }

        int ngap = 0;
        for ( AlignPosList::iterator it = nlist.begin(); it != nlist.end(); ++it ) {
            if ( int(inits[i]+qlen-1) < it->ref_pos ) continue;
            if ( inits[i] >= it->ref_pos ) continue;
            
            if ( param.verbose ) {
                if ( type == DELETION ) {
                    std::cout << "New gap from deletion:";
                    std::cout << reads[i] << "\tinit:" << inits[i] << "\tseq:" << it->ref_pos-inits[i] << "\tref:" << it->ref_pos << "\n";
                }
                else  {
                    std::cout << "New gap from insertion:";
                    std::cout << reads[i] << "\tinit:" << inits[i] << "\tseq:" << it->ref_pos-inits[i]-ngap << "\tref:" << it->ref_pos << "\n";
                }
            }

            // Deletion
            // Sbjct: ABCDEFGHIJKL
            // Query:  BCD---HIJKL
            //  read:   CD---HIJ 
            if ( type == DELETION )
                dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i], it->ref_pos) );

            // Insertion
            // Query:  BCDEFGHIJKL
            // Sbjct: ABCD---HIJKL
            //  read:   CD---HIJ 
            if ( type == INSERTION )
                dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i]-ngap, it->ref_pos) );
            
            ngap++;
            qlen++;
        }    
    }

    /* Is it necessary to put in order??? */
    /* If so, sort them */
    for ( size_t j = 0; j < cdel; j++ ) {
        if ( param.verbose ) std::cout << "Old del:" << dels[j].read << "\t" << dels[j].spos << "\t" << dels[j].rpos << "\n";
        dlist.push_front(dels[j]);
    }

    for ( size_t j = 0; j < cins; j++ ) {
        if ( param.verbose ) std::cout << "Old ins:" << inss[j].read << "\t" << inss[j].spos << "\t" << inss[j].rpos << "\n";
        ilist.push_front(inss[j]);
    }
    
    if ( cdel > 0 ) {
        delete[] dels; cdel = 0;
    }

    if ( cins > 0 ) {
        delete[] inss; cins = 0;
    }

    __listToArray(dlist, DELETION);
    __listToArray(ilist, INSERTION);

    return bad_index;
}

// void SpaPath::insertGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type, Param &param )
// {
//     std::tr1::unordered_map<ReadId, std::list<int> > delmap, insmap;
//     makeMismatchMap(delmap, DELETION);
//     makeMismatchMap(insmap, INSERTION);

//     std::list<Mismatch> dlist, ilist;

//     for ( size_t i = range.first; i < range.second; i++ ) {
        
//         std::string query = bstrs[reads[i]].toString();
//         unsigned     qlen = query.size();
        

//         // Read align start position is greater than the reference position.
//         // Then, increment start position.
//         //     .   
//         // ABCDEF--GHI
//         //         GHI
//         // start = 6 -> 8 
       
//         //        if ( type == DELETION) { 

//         if ( type == INSERTION ) nlist.reverse();
//         for ( AlignPosList::iterator it = nlist.begin(); it != nlist.end(); ++it ) {
//             if ( int(inits[i]+qlen-1) < it->ref_pos ) 
//                 continue;
//             if ( inits[i] >= it->ref_pos ) {
//                 inits[i]++; continue;
//             }
//             //if ( insmap.find(reads[i]) != insmap.end() ) {
//             if ( insmap.find(reads[i]) != insmap.end() ||
//                  delmap.find(reads[i]) != delmap.end() ) {
//                 if ( insmap.find(reads[i]) != insmap.end() ) 
//                     std::cout << "Indel conflict:" << reads[i] << "\n";
//                 if ( delmap.find(reads[i]) != delmap.end() ) 
//                     std::cout << "Deldel check:" << reads[i] << "\n";

//                 dropMismatchRead(i);
                
// //                 int s = inits[i]-qlen; 
// //                 if ( s < 0 ) s = 0;
// //                 int e = inits[i]+2*qlen;
// //                 if ( e >= strlen(consensus) ) e =  strlen(consensus)-1;
//                 int s = 0;
//                 int e = strlen(consensus)-1;
//                 AlignSummary summary = __realignRead(query, s, e);
//                         //__updateIndel(i, inits[i], summary, ilist, dlist);
//                         __updateIndel(i, inits[i], summary, ilist, dlist);
                
//                 //                     appendIndels( summary.ilist, INSERTION );
//                 //                     appendIndels( summary.dlist, DELETION );
//                 //continue;
//                 break;
//             }

// //             for ( int j = 0; j < cdel; j++ ) {
// //                 if ( dels[j].read != reads[i] ) continue;
// //                 std::cout << "Deldel:" << reads[i]<< "\t" << dels[j].spos << "\t" << dels[j].rpos << "\n";
                
// //                 if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
// //             }

// //             for ( int j = 0; j < cins; j++ ) {
// //                 if ( inss[j].read != reads[i] ) continue;
// //                 if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
// //             }

//             std::cout << "New gap:" << reads[i] << "\tinit:" << inits[i] << "\tseq:" << it->ref_pos-inits[i] << "\tref:" << it->ref_pos << "\n";
//             dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i], it->ref_pos) );
//             qlen++;
//         }
//     }
// //     else {
// //         std::vector<int> oinits = std::vector<int>(inits, inits+nread);
// //             bool skip = false;
// //             for ( AlignPosList::reverse_iterator it = nlist.rbegin(); it != nlist.rend(); ++it ) {
// //                 if ( inits[i] > it->ref_pos ) {
// //                     inits[i]++; 
// //                     continue;
// //                 }
// //                 if ( insmap.find(reads[i]) != insmap.end() ) {
// //                     std::cout << "Indel conflict\n";
// //                     dropMismatchRead(i);
// //                     AlignSummary summary = __realignRead(query, inits[i]);
// //                     __updateIndel(i, inits[i], summary, ilist, dlist);
// // //                     appendIndels( summary.ilist, INSERTION );
// // //                     appendIndels( summary.dlist, DELETION );
// //                     skip = true;
// //                     break;
// //                 }

// //                 for ( int j = 0; j < cdel; j++ ) {
// //                     if ( dels[j].read != reads[i] ) continue;
// //                     if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
// //                 }
// //                 for ( int j = 0; j < cins; j++ ) {
// //                     if ( inss[j].read != reads[i] ) continue;
// //                     if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
// //                 }
// //             }

// //             if ( skip ) continue;
            
// //             for ( AlignPosList::reverse_iterator it = nlist.rbegin(); it != nlist.rend(); ++it ) {
// //                 if ( oinits[i] > it->ref_pos ) continue;
// //                 if ( int(oinits[i]+qlen-1) < it->ref_pos ) continue;
// //                 std::cout << "New gap to self:" << reads[i] << "\told-init:" << oinits[i] << "\tnew-init:" << inits[i] << "\tref:" << it->ref_pos << "\n";
// //                 dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i], it->ref_pos) );
// //                 //qlen++;
// //             }

// //         }
// //     }

//     /* Is it necessary to put in order??? */
//     /* If so, sort them */
//     for ( int j = 0; j < cdel; j++ ) 
//         dlist.push_front(dels[j]);

//     for ( int j = 0; j < cins; j++ ) 
//         ilist.push_front(inss[j]);
    
//     if ( cdel > 0 ) {
//         delete[] dels; 
//         cdel = 0;
//     }

//     if ( cins > 0 ) {
//         delete[] inss; 
//         cins = 0;
//     }

//     __listToArray(dlist, DELETION);
//     __listToArray(ilist, INSERTION);

// }


// void SpaPath::insertGapsToReads( AlignPosList &nlist, BitString *bstrs, int type )
// {
//     std::list<Mismatch> dlist;
// //     if ( cdel > 0 ) {
// //         dlist = std::list<Mismatch>( dels, dels+cdel );
// //         delete[] dels; 
// //         cdel = 0;
// //     }

//     //for ( AlignPosList::reverse_iterator it = nlist.rbegin(); it != nlist.rend(); ++it ) {
//     //for ( AlignPosList::iterator it = nlist.begin(); it != nlist.end(); ++it ) {
//     //std::cout << "Gap pos:" << it->ref_pos << "\n";

//     for ( size_t i = 0; i < nread; i++ ) {
        
//         std::string query = bstrs[reads[i]].toString();
//         unsigned     qlen = query.size();
        

//         // Read align start position is greater than the reference position.
//         // Then, increment start position.
//         //     .   
//         // ABCDEF--GHI
//         //         GHI
//         // start = 6 -> 8 
       
//         if ( type == DELETION) { 
//             for ( AlignPosList::iterator it = nlist.begin(); it != nlist.end(); ++it ) {
//                 //for ( size_t i = 0; i < nread; i++ ) {
//                 if ( inits[i] >= it->ref_pos ) {
//                     inits[i]++; 
//                     continue;
//                 }


//                 for ( int j = 0; j < cdel; j++ ) {
//                     if ( dels[j].read != reads[i] ) continue;
//                     if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
//                 }

//                 for ( int j = 0; j < cins; j++ ) {
//                     if ( inss[j].read != reads[i] ) continue;
//                     if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
//                 }
//                 //             std::string query = bstrs[reads[i]].toString();
//                 //             unsigned     qlen = query.size();
//                 if ( int(inits[i]+qlen-1) < it->ref_pos ) {
//                     continue;
//                 }
                
//                 std::cout << "New gap to other:" << reads[i] << "\tinit:" << inits[i] << "\tseq:" << it->ref_pos-inits[i] << "\tref:" << it->ref_pos << "\n";
//                 dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i], it->ref_pos) );
//                 qlen++;
//             }
//         }
//         else {
//             std::vector<int> oinits = std::vector<int>(inits, inits+nread);
//             for ( AlignPosList::reverse_iterator it = nlist.rbegin(); it != nlist.rend(); ++it ) {
//                 //for ( size_t i = 0; i < nread; i++ ) {
//                 if ( inits[i] > it->ref_pos ) {
//                     inits[i]++; 
//                     continue;
//                 }

//                 for ( int j = 0; j < cdel; j++ ) {
//                     if ( dels[j].read != reads[i] ) continue;
//                     if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
//                 }
//                 for ( int j = 0; j < cins; j++ ) {
//                     if ( inss[j].read != reads[i] ) continue;
//                     if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
//                 }

// //                 //             std::string query = bstrs[reads[i]].toString();
// //                 //             unsigned     qlen = query.size();
// //                 if ( int(inits[i]+qlen-1) < it->ref_pos ) {
// //                     continue;
// //                 }
                
//             }

//             for ( AlignPosList::reverse_iterator it = nlist.rbegin(); it != nlist.rend(); ++it ) {
//                 if ( oinits[i] > it->ref_pos ) continue;
//                 if ( int(oinits[i]+qlen-1) < it->ref_pos ) continue;
//                 std::cout << "New gap to self:" << reads[i] << "\told-init:" << oinits[i] << "\tnew-init:" << inits[i] << "\tref:" << it->ref_pos << "\n";
//                 dlist.push_back( Mismatch(reads[i], it->ref_pos-inits[i], it->ref_pos) );
//                 //qlen++;
//             }

//         }
//     }


//     /* Is it necessary to put in order??? */
//     /* If so, sort them */
//     for ( int j = 0; j < cdel; j++ ) 
//         dlist.push_front(dels[j]);
    
//     if ( cdel > 0 ) {
//         delete[] dels; 
//         cdel = 0;
//     }

//     __listToArray(dlist, DELETION);
// }

// Join two objects
// After alignment
void SpaPath::join( SpaPath *oaln, int start, int lgap, int tgap, AlignPosList &ninss, AlignPosList &ndels, int kmer_size, InvertedIndex &iindex, BitString *bstrs, Param &param )
{
    if ( lgap ) prepend(oaln, lgap, kmer_size);
    else if ( tgap ) append(oaln, tgap, kmer_size);

    if ( lgap ) adjustPositions(lgap);
    else oaln->adjustPositions(start);

    if ( param.verbose ) {
        std::cout << "#insertions:" << cins << "\n";
        std::cout << "#deletions: " << cdel << "\n";
    }

    int cnread = nread;
    //int onread = oaln->getReadCount();

    mergeReads(oaln);
    mergeIndels(oaln);
    
    if ( param.verbose ) {
        std::cout << "#insertions:" << cins << "\n";
        std::cout << "#deletions: " << cdel << "\n";
    }

    /* no indels between two paths */
    if ( ninss.size() == 0 && ndels.size() == 0 ) {
        if (param.verbose) std::cout << "No indels - Skip indel adjustment\n";
        //mergeReads(oaln);
        //mergeIndels(oaln);
        return;
    }

    /* Both insertion and deletion 
       How ??
    */
    if ( ndels.size() && ninss.size() ) {
        if ( param.verbose ) std::cout << "Both insertion/deletion\n";
        //return;
    }


    /* Insertion
       ABCD----IJKLM
       ABCDEFGHIJK

       Update consensus with gaps.
       Update start of sbjct reads.
       Update sbjct position of gaps.
       For each sbjct read, insert gaps if no carry-over indels. 
       Otherwise, re-align to gapppy consensus.
       //For each query read with indel, realign to consensus (gappy).
       Leave intact of indels from query reads.
     */

    makeGappyConsensus( ninss, param );

    std::list<int> bad_index;
    if ( ninss.size() ) {
        bad_index = updateGaps( ninss, bstrs, IntPair(0, cnread), INSERTION, param );
        dropReads(bad_index);
        //return;
    }

    /* Deletion 
       ABCDEFGHIJKLM
       ABCD----IJK
       
       Update start of reads from query.
       
       Leave intact of indels from sbject reads.
       For each query read, insert gaps if no carry-over indels. 
       Otherwise, re-align to consensus.
    */
    
    if ( ndels.size() ) {        
        bad_index = updateGaps( ndels, bstrs, IntPair(cnread, nread), DELETION, param );
        dropReads(bad_index);
        //return;
    }
    

//     for ( AlignPosList::reverse_iterator it = ndels.rbegin(); it != ndels.rend(); ++it ) {
//         for ( size_t i = cnread; i < nread; i++ ) {
//             for ( int j = 0; j < cins; j++ ) {
//                 if ( inss[j].read != reads[i] ) continue;
//                 if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
//             }
//             for ( int j = 0; j < cdel; j++ ) {
//                 if ( dels[j].read != reads[i] ) continue;
//                 if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
//             }
//         }
//     }
    


//     std::string nstr = biostr::getSequenceString( kmers, nkmer, kmer_size );
//     if ( param.verbose ) std::cout << "Consensus:" << nstr << "\n";
//     for ( AlignPosList::reverse_iterator it = ninss.rbegin(); it != ninss.rend(); ++it ) {
//         assert(it->ref_pos < nstr.size());
//         nstr.insert( it->ref_pos, "-" );  
//         if ( param.verbose ) std::cout << it->ref_pos << "\n";  

//         for ( size_t i = 0; i < cnread; i++ ) {
//             for ( int j = 0; j < cins; j++ ) {
//                 if ( inss[j].read != reads[i] ) continue;
//                 if ( inss[j].rpos >= it->ref_pos ) inss[j].rpos++;
//             }
//             for ( int j = 0; j < cdel; j++ ) {
//                 if ( dels[j].read != reads[i] ) continue;
//                 if ( dels[j].rpos >= it->ref_pos ) dels[j].rpos++;
//             }
//         }
//     }

//     if ( param.verbose) std::cout << "New Consensus:" << nstr << "\n";
//     KmerArray nkmers = biostr::getKmers(nstr, kmer_size);
//     updateKmers( &nkmers[0], nkmers.size() );
//     //updateConsensusSequence( nstr.c_str() );

//     for ( size_t i = 0; i < nkmers.size(); i++ ) 
//         std::cout << alpha::IntegerToAminoAcid(nkmers[i], kmer_size) << " ";
//     std::cout << "\n";

//     if (param.verbose) std::cout << "Deletion adjustment\n";
//     if ( ninss.size() ) std::cout << "Insert gaps from insertion\n";
//     for ( AlignPosList::iterator it = ninss.begin(); it != ninss.end(); ++it ) 
//         std::cout << it->seq_pos << ":" << it->ref_pos << "\n";
//     std::cout << "\n";
//     if ( ninss.size() ) insertGaps( ninss, bstrs, IntPair(0, cnread), INSERTION );
//     if ( ninss.size() ) verifyGaps( ninss, bstrs, IntPair(cnread, nread), INSERTION );
//     if ( ndels.size() ) std::cout << "Insert gaps from deletion\n";
//     for ( AlignPosList::iterator it = ndels.begin(); it != ndels.end(); ++it ) 
//         std::cout << it->seq_pos << ":" << it->ref_pos << "\n";
//     std::cout << "\n";
//     if ( ndels.size() ) insertGaps( ndels, bstrs, IntPair(cnread, nread), DELETION );
//     if ( ndels.size() ) verifyGaps( ndels, bstrs, IntPair(0, cnread), DELETION );
// //     if ( ninss.size() ) insertGapsToReads(ninss, bstrs, INSERTION);
// //     if ( ndels.size() ) oaln->insertGapsToReads(ndels, bstrs, DELETION);
    
// //     mergeReads(oaln);
// //     mergeIndels(oaln);

//     //if ( cins > 0  ) delete[] inss;
//     //if ( cdel > 0  ) delete[] dels;
//     //cins = cdel = 0;
//     //inss = dels = NULL;

// //     if (param.verbose) std::cout << "Align reads to path\n";
// //     align(bstrs, param);
// //     validate(bstrs, param);
}


void SpaPath::join( SpaPath *oaln, std::string &mid, int start, int lgap, int tgap, AlignPosList &ninss, AlignPosList &ndels, int kmer_size,  InvertedIndex &iindex, BitString *bstrs, Param &param )
{
    if ( param.verbose ) std::cout << "Middle:" << mid << "\n";
    append(mid, kmer_size);
    if ( lgap )   prepend(oaln, lgap, kmer_size);
    else if ( tgap )   append(oaln, tgap, kmer_size);

    if ( lgap ) adjustPositions(lgap);
    else oaln->adjustPositions(start);

    if ( param.verbose ) {
        std::cout << "#insertions:" << cins << "\n";
        std::cout << "#deletions: " << cdel << "\n";
    }

    int cnread = nread;
    //int onread = oaln->getReadCount();

    mergeReads(oaln);
    mergeIndels(oaln);
    
    if ( param.verbose ) {
        std::cout << "#insertions:" << cins << "\n";
        std::cout << "#deletions: " << cdel << "\n";
    }

    /* no indels between two paths */
    if ( ninss.size() == 0 && ndels.size() == 0 ) {
        if (param.verbose) std::cout << "No indels - Skip indel adjustment\n";
        //mergeReads(oaln);
        //mergeIndels(oaln);
        return;
    }

    /* Both insertion and deletion 
       How ??
    */
    if ( ndels.size() && ninss.size() ) {
        std::cout << "Both insertion/deletion: check correctedness\n";
        //return;
    }


    /* Insertion
       ABCD----IJKLM
       ABCDEFGHIJK

       Update consensus with gaps.
       Update start of sbjct reads.
       Update sbjct position of gaps.
       For each sbjct read, insert gaps if no carry-over indels. 
       Otherwise, re-align to gapppy consensus.
       //For each query read with indel, realign to consensus (gappy).
       Leave intact of indels from query reads.
     */

    makeGappyConsensus( ninss, param );

    std::list<int> bad_index;
    if ( ninss.size() ) {
        bad_index = updateGaps( ninss, bstrs, IntPair(0, cnread), INSERTION, param );
        dropReads(bad_index);
        //return;
    }

    /* Deletion 
       ABCDEFGHIJKLM
       ABCD----IJK
       
       Update start of reads from query.
       
       Leave intact of indels from sbject reads.
       For each query read, insert gaps if no carry-over indels. 
       Otherwise, re-align to consensus (ungappy).
    */
    
    if ( ndels.size() ) {        
        bad_index = updateGaps( ndels, bstrs, IntPair(cnread, nread), DELETION, param );
        dropReads(bad_index);
        //return;
    }

    
//     std::cout << "#insertions:" << cins << "\n";
//     std::cout << "#deletions: " << cdel << "\n";
    
//     int cnread = nread;
//     //int onread = oaln->getReadCount();

//     mergeReads(oaln);
//     mergeIndels(oaln);

//     std::cout << "#insertions:" << cins << "\n";
//     std::cout << "#deletions: " << cdel << "\n";
    

//     /* no indels between two paths */
//     if ( ninss.size() == 0 && ndels.size() == 0 ) {
//         if (param.verbose) std::cout << "No indels - Skip indel adjustment\n";
// //         mergeReads(oaln);
// //         mergeIndels(oaln);
//         return;
//     }

//     std::string nstr = biostr::getSequenceString( kmers, nkmer, kmer_size );
//     for ( AlignPosList::reverse_iterator it = ninss.rbegin(); it != ninss.rend(); ++it ) {
//         assert(it->ref_pos < nstr.size());
//         nstr.insert( it->ref_pos, "-" );    
//     }
//     if ( param.verbose) std::cout << "New Consensus:" << nstr << "\n";
//     KmerArray nkmers = biostr::getKmers(nstr, kmer_size);
//     updateKmers( &nkmers[0], nkmers.size() );
//     //updateConsensusSequence( nstr.c_str() );

//     if (param.verbose) std::cout << "Deletion adjustment\n";
//     if ( ninss.size() ) std::cout << "Insert gaps from insertion\n";
//     for ( AlignPosList::iterator it = ninss.begin(); it != ninss.end(); ++it ) 
//         std::cout << it->seq_pos << ":" << it->ref_pos << "\n";
//     std::cout << "\n";
//     if ( ninss.size() ) insertGaps( ninss, bstrs, IntPair(0, cnread), INSERTION );
//     if ( ndels.size() ) std::cout << "Insert gaps from deletion\n";
//     for ( AlignPosList::iterator it = ndels.begin(); it != ndels.end(); ++it ) 
//         std::cout << it->seq_pos << ":" << it->ref_pos << "\n";
//     std::cout << "\n";
//     if ( ndels.size() ) insertGaps( ndels, bstrs, IntPair(cnread, nread), DELETION );
// //     if ( ninss.size() ) insertGapsToReads(ninss, bstrs, INSERTION);
// //     if ( ndels.size() ) oaln->insertGapsToReads(ndels, bstrs, DELETION);
    
// //     mergeReads(oaln);
// //     mergeIndels(oaln);

    
// //     //if ( cins > 0  ) delete[] inss;
// //     //if ( cdel > 0  ) delete[] dels;
// //     //cins = cdel = 0;
// //     //inss = dels = NULL;

// //     if (param.verbose) std::cout << "Align reads to path\n";
// //     align(bstrs, param);
// //     validate(bstrs, param);
}

// void SpaPath::join( SpaPath *oaln, KmerId *mkmers, int mnkmer, int start, int lgap, int tgap, AlignPosList &ilist, AlignPosList &dlist, int kmer_size )
// {
//     if ( mnkmer ) append(mkmers, mnkmer);
//     if ( lgap )   prepend(oaln, lgap);
//     else if ( tgap )   append(oaln, tgap, kmer_size);

//     if ( lgap ) adjustPositions(lgap);
//     else oaln->adjustPositions(start);

//     mergeReads(oaln);

    
//     // indel handling of two consensus sequences
//     // if insertion, put '-' in all reads in ref
//     // if deletion, put '-' in all reads in query
//     //mergeIndel(ilist, dlist);
// }

std::pair<int, PathId> SpaPath::getLeadingGapCount( PathAlignPairList &path_aligns )
{
    int lgap = 0;
    PathId pid = NOT_PATH;
    for ( PathAlignPairList::iterator jt = path_aligns.begin(); jt != path_aligns.end(); ++jt ) {
        PathId pi = jt->first;
        int lg = jt->second.lgap.first;
        if ( lg > lgap ) {
            pid = pi;
            lgap = lg;
        }
    }
    return std::pair<int, PathId>(lgap, pid);
}

std::pair<int, PathId> SpaPath::getTrailingGapCount( PathAlignPairList &path_aligns )
{
    int tgap = 0;
    PathId pid = NOT_PATH;
    for ( PathAlignPairList::iterator jt = path_aligns.begin(); jt != path_aligns.end(); ++jt ) {
        PathId pi = jt->first;
        int tg = jt->second.egap.first;
        if ( tg > tgap ) {
            pid = pi;
            tgap = tg;
        }
    }
    return std::pair<int, PathId>(tgap, pid);
}

void SpaPath::mergeCluster( PathAlignPairList &path_aligns,
                            std::tr1::unordered_map<PathId, SpaPath*> &path2aln_map,
                            BitString *bstrs,
                            Param &param )
{
    std::pair<int, PathId> lgap = getLeadingGapCount( path_aligns );
    std::pair<int, PathId> tgap = getTrailingGapCount( path_aligns );
    
    if ( param.verbose ) std::cout << "lgap:" << lgap.first << "\tpath:" << lgap.second << "\n";
    if ( param.verbose ) std::cout << "tgap:" << tgap.first << "\tpath:" << tgap.second << "\n";

    //if ( param.verbose ) std::cout << consensus << "\n";
    if ( lgap.first )   prepend(path2aln_map[lgap.second], lgap.first, param.kmer_size);
    //if ( param.verbose ) std::cout << consensus << "\n";
    if ( tgap.first )   append(path2aln_map[tgap.second], tgap.first, param.kmer_size);
    //if ( param.verbose ) std::cout << consensus << "\n";

    /* adjust read start/indel for this path */
    if ( lgap.first ) adjustPositions(lgap.first);
      

    PathAlignPairList ins_pairs, del_pairs;
    for ( PathAlignPairList::iterator it = path_aligns.begin(); it != path_aligns.end(); ++it ) {
        PathId query_pid = it->first;
        AlignSummary summary = it->second;
        if ( summary.ilist.size() ) {
            ins_pairs.push_back( *it );
            continue;
        }
        else if ( summary.dlist.size() ) {
            del_pairs.push_back( *it );
            continue;
        }

        int start = 0;
        if ( summary.lgap.first > 0 ) 
            start = lgap.first - summary.lgap.first;
        else start = summary.range.first + lgap.first;

        if ( param.verbose ) std::cout << "query:" << query_pid 
                                       << "\tlgap:" << summary.lgap.first << "," << summary.lgap.second
                                       << "\tegap:" << summary.egap.first << "," << summary.egap.second
                                       << "\trange:" << summary.range.first << "," << summary.range.second << "\t"
                                       << "\tstart:" << start << "\n";

        path2aln_map[query_pid]->adjustPositions(start);
        mergeReads(path2aln_map[query_pid]);
        mergeIndels(path2aln_map[query_pid]);
    }

    if ( param.verbose ) std::cout << "Alignment w/ insertions:" << ins_pairs.size() << "\n";
    AlignPosList ninss;
    for ( PathAlignPairList::iterator it = ins_pairs.begin(); it != ins_pairs.end(); ++it ) {
        PathId query_pid = it->first;
        AlignSummary summary = it->second;
        int start = 0;
        if ( summary.lgap.first > 0 ) 
            start = lgap.first - summary.lgap.first;
        else start = summary.range.first + lgap.first;

        if ( param.verbose ) std::cout << "query:" << query_pid 
                                       << "\tlgap:" << summary.lgap.first << "," << summary.lgap.second
                                       << "\tegap:" << summary.egap.first << "," << summary.egap.second
                                       << "\trange:" << summary.range.first << "," << summary.range.second << "\t"
                                       << "\tstart:" << start << "\n";

        path2aln_map[query_pid]->adjustPositions(start);
        mergeReads(path2aln_map[query_pid]);
        mergeIndels(path2aln_map[query_pid]);

        for ( AlignPosList::iterator jt = summary.ilist.begin(); jt != summary.ilist.end(); ++jt )
            ninss.push_back(*jt);
    }

    size_t cnread = nread;

    std::list<int> bad_index;
    if ( ninss.size() ) {
        makeGappyConsensus( ninss, param );
        bad_index = updateGaps( ninss, bstrs, IntPair(0, cnread), INSERTION, param );
        dropReads(bad_index);
    }

    if ( param.verbose ) std::cout << "Alignment w/ deletion:" << del_pairs.size() << "\n";
    AlignPosList ndels;
    for ( PathAlignPairList::iterator it = del_pairs.begin(); it != del_pairs.end(); ++it ) {
        PathId query_pid = it->first;
        AlignSummary summary = it->second;
        int start = 0;
        if ( summary.lgap.first > 0 ) 
            start = lgap.first - summary.lgap.first;
        else start = summary.range.first + lgap.first;

        if ( param.verbose ) std::cout << "query:" << query_pid 
                                       << "\tlgap:" << summary.lgap.first << "," << summary.lgap.second
                                       << "\tegap:" << summary.egap.first << "," << summary.egap.second
                                       << "\trange:" << summary.range.first << "," << summary.range.second << "\t"
                                       << "\tstart:" << start << "\n";

        path2aln_map[query_pid]->adjustPositions(start);
        mergeReads(path2aln_map[query_pid]);
        mergeIndels(path2aln_map[query_pid]);

        for ( AlignPosList::iterator jt = summary.ilist.begin(); jt != summary.ilist.end(); ++jt )
            ndels.push_back(*jt);
    }

    if ( ndels.size() ) {        
        bad_index = updateGaps( ndels, bstrs, IntPair(cnread, nread), DELETION, param );
        dropReads(bad_index);
    }

    
//         for ( AlignPosList::iterator jt = summary.ilist.begin(); jt != summary.ilist.end(); ++jt )
//             ninss.push_back(*jt);
//         for ( AlignPosList::iterator jt = summary.dlist.begin(); jt != summary.dlist.end(); ++jt )
//             ndels.push_back(*jt);
    
//         if ( param.verbose ) {
//             std::cout << "#insertions:" << cins << "\n";
//             std::cout << "#deletions: " << cdel << "\n";
//         }
        
//         /* no indels between two paths */
//         if ( ninss.size() == 0 && ndels.size() == 0 ) {
//             if (param.verbose) std::cout << "No indels - Skip indel adjustment\n";
//             continue;
//         }
        
//         if ( ndels.size() && ninss.size() ) {
//             std::cout << "Both insertion/deletion: check correctedness\n";
//         }
        
        
//         /* Insertion
//            ABCD----IJKLM
//            ABCDEFGHIJK
           
//            Update consensus with gaps.
//            Update start of sbjct reads.
//            Update sbjct position of gaps.
//            For each sbjct read, insert gaps if no carry-over indels. 
//            Otherwise, re-align to gapppy consensus.
//            //For each query read with indel, realign to consensus (gappy).
//            Leave intact of indels from query reads.
//         */
        


//     std::list<int> bad_index;
//     if ( ninss.size() ) {
//         bad_index = updateGaps( ninss, bstrs, IntPair(0, nread), INSERTION, param );
//         dropReads(bad_index);
//         //return;
//     }

//     /* Deletion 
//        ABCDEFGHIJKLM
//        ABCD----IJK
       
//        Update start of reads from query.
       
//        Leave intact of indels from sbject reads.
//        For each query read, insert gaps if no carry-over indels. 
//        Otherwise, re-align to consensus (ungappy).
//     */
    
//     if ( ndels.size() ) {        
//         bad_index = updateGaps( ndels, bstrs, IntPair(0, nread), DELETION, param );
//         dropReads(bad_index);
//         //return;
//     }
}


// // 
// void SpaPath::latch( SpaPath *raln, int start )
// {
//     unsigned rnkmer = raln->getKmerCount();
//     unsigned rkmers = raln->getKmers();

//     size_t old_nkmer  = nkmer;
//     KmerId *old_kmers = kmers;

//     int diff = nkmer - start;
//     if ( diff >= 0 ) {
//         nkmer += (rnkmer-diff);            
//     } else {
//         nread += (-1*diff+rnkmer);
//     }

//     kmers = new KmerId[nkmer];
//     for ( size_t i = 0; i < old_nkmer; i++ )
//         kmers[i] = old_kmers[i];
//     if ( diff >= 0 ) {
//         nkmer += (rnkmer-diff);            
//     } else {
//         nread += (-1*diff+rnkmer);
//     }
    
// }

void SpaPath::addRead( ReadId rid, int spos )
{
    size_t old_nreads = nread;
    ReadId *old_reads = reads;
    int *old_inits = inits;

    nread += 1;
    reads = new ReadId[nread];
    inits = new int[nread];
    for ( size_t i = 0; i < old_nreads; i++ ){
        reads[i] = old_reads[i];
        inits[i] = old_inits[i];
    }
    reads[old_nreads] = rid;
    inits[old_nreads] = spos;
    
    delete[] old_reads;
    delete[] old_inits;
}

void SpaPath::addReads( ReadId *rids, int *sposs, size_t nrid )
{
    size_t old_nreads = nread;
    ReadId *old_reads = reads;
    int *old_inits = inits;

    nread += nrid;
    reads = new ReadId[nread];
    inits = new int[nread];
    for ( size_t i = 0; i < old_nreads; i++ ){
        reads[i] = old_reads[i];
        inits[i] = old_inits[i];
    }
    for ( size_t i = 0; i < nrid; i++ ){
        reads[i+old_nreads] = rids[i];
        inits[i+old_nreads] = sposs[i];
    }
    delete[] old_reads;
    delete[] old_inits;
}

void SpaPath::joinReads( ReadId *rids, int *poss, int rsize, BitString *bstrs, int k, int direction, bool verbose )
{
    direction == 0 ? __joinReadsLeft( rids, poss, rsize, bstrs, k, verbose ) : __joinReadsRight( rids, poss, rsize, bstrs, k, verbose );
}

void SpaPath::__joinReadsLeft( ReadId *rids, int *sposs, int rsize, BitString *bstrs, int k, bool verbose)
{
    int max = 0;
    for ( int i = 0; i < rsize; i++ ) 
        if ( sposs[i] > max ) max = sposs[i];

    if (verbose) std::cout << "Max:" << max << "\n";

    std::string ostr = biostr::getSequenceString(kmers, nkmer,k);
    std::string nstr = std::string(max, 'X');
    nstr = nstr + ostr;
    if (verbose) std::cout << "New:" << nstr << "\n";
    KmerArray nkids = biostr::getKmers(nstr, k);

    KmerId *okmers = kmers;
    unsigned osize = nkmer;
    nkmer = nkids.size();
    kmers = new KmerId[nkmer];

    for ( unsigned i = 0; i < nkmer; i++ ) kmers[i] = nkids[i];
    delete[] okmers;

    adjustPositions(max);
    
    ReadId *oreads = reads;
    osize = nread;
    nread += rsize;
    reads = new ReadId[nread];

    for ( unsigned i = 0; i < osize; i++ ) reads[i] = oreads[i];
    for ( int i = 0; i < rsize; i++ ) reads[i+osize] = rids[i];
    delete[] oreads;

    int *oinit = inits;
    inits = new int[nread];

    for ( unsigned i = 0; i < osize; i++ ) inits[i] = oinit[i];
    for ( int i = 0; i < rsize; i++ ) inits[osize+i] = max-sposs[i];    
    delete[] oinit;
}

void SpaPath::__joinReadsRight( ReadId *rids, int *rposs, int rsize, BitString *bstrs, int k, bool verbose )
{
    int slen = nkmer+k-1;

    int max = 0;
    for ( int i = 0; i < rsize; i++ ) { 
        int rlen = bstrs[rids[i]].toString().length();
        int npos = rlen - (slen-rposs[i]);
        if ( npos > max ) max = npos;
    }
    
    if (verbose) std::cout << "Max:" << max << "\n";

    std::string ostr = biostr::getSequenceString(kmers, nkmer,k);
    std::string nstr = std::string(max, 'X');
    nstr = ostr + nstr;

    if ( verbose) std::cout << "New:" << nstr << "\n";
    KmerArray nkids = biostr::getKmers(nstr, k);

    KmerId *okmers = kmers;
    unsigned osize = nkmer;
    nkmer = nkids.size();
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ ) kmers[i] = nkids[i];
    delete[] okmers;

    //adjustPositions(max);
    
    ReadId *oreads = reads;
    osize = nread;
    nread += rsize;
    reads = new ReadId[nread];

    for ( unsigned i = 0; i < osize; i++ ) reads[i] = oreads[i];
    for ( int i = 0; i < rsize; i++ ) reads[i+osize] = rids[i];
    delete[] oreads;

    int *oinit = inits;
    inits = new int[nread];

    for ( unsigned i = 0; i < osize; i++ ) inits[i] = oinit[i];
    for ( int i = 0; i < rsize; i++ ) inits[osize+i] = rposs[i];    
    delete[] oinit;
}

// void SpaPath::updateConsensus( std::string &str, BitString *bstrs, double merge_score, int k, bool verbose )
// {
//     KmerArray kids = biostr::getKmers(str, k);
// //     if ( nkmer != kids.size() ) {
// //         std::cout << "New kmer size is not identical\n"; exit(1);
// //     }
//     bool same = true;
//     if ( nkmer != kids.size() ) same = false;

//     if ( same ) {
//         for ( size_t i = 0; i < nkmer; i++ ) kmers[i] = kids[i];
//     } else {
//         KmerId *okmers = kmers;
//         unsigned osize = nkmer;
//         nkmer = kids.size();
//         kmers = new KmerId[nkmer];
//         for ( size_t i = 0; i < nkmer; i++ ) kmers[i] = kids[i];
//         delete[] okmers;
        
//         std::string nstr = biostr::stripGap( biostr::getSequenceString(kmers, nkmer, k) );
//         std::string ostr = biostr::stripGap( biostr::getSequenceString(okmers, osize, k) );
        
//         __updateReadPos(nstr, bstrs, merge_score, verbose);
//     }
// }

// void SpaPath::__updateReadPos( std::string &consensus, BitString *bstrs, double merge_score, bool verbose )
// {
//     uList mis;
//     std::list<Mismatch>ilist, dlist; // temprary lists;    
//     for ( unsigned i = 0; i < nread; i++ ) {

//         std::string seq = bstrs[reads[i]].toString();
//         unsigned seq_len = seq.size();
//         unsigned ref_len = consensus.size();

//         if (verbose) {
//             std::cout << "\t" << i << "\t" << reads[i] << "\t";
//             std::cout << "seq len:" << seq_len << "\tref len:" << ref_len << "\n";
//             std::cout << i << "\tseq:" << seq << "\n";
//         }
        
//         int s = inits[i];
//         int l = seq_len;        
//         if ( verbose ) std::cout << "s:" << s << "\tl:" << l << "\n";
//         if ( s < 0 ) {
//             l += s;
//             s = 0;
//         }
//         if ( verbose ) std::cout << "s:" << s << "\tl:" << l << "\n";


//         //int len = seq_len;
//         if ( s+l > (int)ref_len ) l -= (s+l - ref_len);
//         if ( verbose ) std::cout << "s:" << s << "\tl:" << l << "\n";


//         assert( s+l <= (int)consensus.size() );
//         std::string ref = consensus.substr(s, l);
//         if ( verbose ) std::cout << i << "\tref:" << ref << "\n";

//         int pos = scoring::countPositive(seq, ref, BLOSUM62);
//         assert(seq_len > 0);
//         double rate = (double)pos/seq_len;
//         if ( verbose ) std::cout << "Positive score:" << rate << "\n";
//         if ( rate >= merge_score ) {
//             continue;
//         }

//         GlobalAlignPair paln = GlobalAlignPair(ref, seq);
//         AlignSummary summary = paln.getSummary();

//         size_t aln_beg = summary.range.first;
//         size_t aln_end = summary.range.second;
//         size_t aln_len = aln_end - aln_beg + 1;
                    
//         if ( verbose ) {
//             std::cout << paln.getAlignment();
            
//             std::cout << "\t#ins:" << summary.ins.size() << "\t";
//             std::cout << "\t#del:" << summary.del.size() << "\t";
//             std::cout << "\t#mat:" << summary.match << "\t";
//             std::cout << "\t#pos:" << summary.positive << "\t";
//             std::cout << "\t#mms:" << summary.mismatch << "\n";
//             std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
//             std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
//         }

//         int rpos = inits[i];

        
//         if ( summary.lgap.first == 0 && summary.lgap.second > 0 ) {
//             if ( verbose ) std::cout << "\tstart:" << inits[i];
//             if ( verbose ) std::cout << " -> ";
//             int diff = (summary.lgap.second-summary.lgap.first);
//             inits[i] = inits[i] + diff;
//             if ( summary.range.first != diff ) inits[i] += summary.range.first;
//             if ( verbose ) std::cout << inits[i];
//             if ( verbose ) std::cout << "\n";
//         } else if ( summary.lgap.first > 0 && summary.lgap.second == 0 ) {
//             if ( verbose ) std::cout << "\tstart:" << inits[i];
//             if ( verbose ) std::cout << " -> ";
//             int diff =(summary.lgap.first-summary.lgap.second);
//             inits[i] = inits[i] - diff;
//             if ( summary.range.first != diff ) inits[i] += summary.range.first;
//             if ( verbose ) std::cout << inits[i];
//             if ( verbose ) std::cout << "\n";
//         }
        

//         if ( aln_len < unsigned(merge_score * seq_len) ) {
//             if ( verbose ) std::cout << "\talignment - short\n";
//             mis.push_back(i);
//             continue;
//         }
//         else if ( aln_len * merge_score  >  seq_len ) {
//             if ( verbose ) std::cout << "\talignment - long\n";
//             mis.push_back(i);
//             continue;
//         }
        
//         //if ( summary.positive < aln_len * pcut ) {
//         if ( (double)summary.positive <  aln_len * merge_score ) {
//             mis.push_back(i);
//             if ( verbose ) std::cout << "\talignment - weak\n";
//             continue;
//         }

//         if ( verbose ) {
//             if ( summary.ilist.size() > 0 ) {
//                 std::cout << "\tInsertions:\t";
//                 for (AlignPosList::iterator it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
//                     //std::cout << "aln:" << it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << it->ref_pos + inits[i] << "\t";
//                     std::cout << "aln:" << it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << rpos + it->ref_pos << "\t";
//                 }
//             } 
//             if ( summary.dlist.size() > 0 ) {
//                 std::cout << "\tDeletions:\t";
//                 for (AlignPosList::iterator it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
//                     //std::cout << "aln:" <<it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << it->ref_pos + inits[i] << "\t";
//                     std::cout << "aln:" <<it->aln_pos << ", seq:" << it->seq_pos << ", ref:" << it->ref_pos << ", adj:" << rpos + it->ref_pos << "\t";
//                 }
//             }
//             std::cout << "\n";
//         }

//         AlignPosList::iterator it;
//         for ( it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
//             ilist.push_back( Mismatch(reads[i], it->seq_pos, rpos + it->ref_pos) );
//         }
//         for ( it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
//             dlist.push_back( Mismatch(reads[i], it->seq_pos, rpos + it->ref_pos) );
//         }
//     }

//     Mismatch *oins = inss;
//     Mismatch *odel = dels;
//     oins = odel = NULL;
//     if ( cins > 0 ) delete[] oins;
//     if ( cdel > 0 ) delete[] odel;
//     cins = cdel = 0;

//     __listToArray(ilist, INSERTION);
//     __listToArray(dlist, DELETION);
//     __trim(mis, verbose);
//     //delete[] oins;
//     //delete[] odel;

// }


void SpaPath::updateKmers( KmerId *nkids, unsigned nsize )
{
    KmerId *okmers = kmers;
    nkmer = nsize;
    
    kmers = new KmerId[nkmer];
    for ( unsigned i = 0; i < nkmer; i++ )
        kmers[i] = nkids[i];
    
    delete[] okmers;
}


void SpaPath::updateConsensusSequence( const char*nstr )
{
    char *ostr = consensus;
    int olen = strlen(consensus);

    int len = strlen(nstr);
    consensus = new char[len+1];
    strcpy(consensus, nstr);
    if ( olen > 0 ) delete[] ostr;
}

void SpaPath::dropMismatchRead( int index )
{
    std::list<Mismatch>ilist, dlist;
//     if ( cins ) ilist = std::list<Mismatch>( inss, inss+cins );
//     if ( cdel ) dlist = std::list<Mismatch>( dels, dels+cdel );

//     std::list<Mismatch>::iterator it;
//     for ( it = ilist.begin(); it != ilist.end(); ) {
//         if ( reads[index] == (*it).read ) ilist.erase(it++);
//         else ++it;
//     } 
//     for ( it = dlist.begin(); it != dlist.end(); ) {
//         if ( reads[index] == (*it).read ) dlist.erase(it++);
//         else ++it;
//     }

    for ( size_t i = 0; i < cins; i++ )
        if ( inss[i].read != reads[index] ) 
            ilist.push_back(inss[i]);

    for ( size_t i = 0; i < cdel; i++ )
        if ( dels[i].read != reads[index] ) 
            dlist.push_back(dels[i]);

    //if ( ilist.size() != cins ) 
    __listToArray(ilist, INSERTION);
    //if ( dlist.size() != cdel )
    __listToArray(dlist, DELETION);
    
}
void SpaPath::dropMismatches( ReadIdArray &rids )
{
    std::list<Mismatch>ilist, dlist;
    ilist = std::list<Mismatch>( inss, inss+cins );
    dlist = std::list<Mismatch>( dels, dels+cdel );

    std::list<Mismatch>::iterator it;
    for ( size_t i = 0; i < rids.size(); i++ ) {
        for ( it = ilist.begin(); it != ilist.end(); ) {
            if ( rids[i] == (*it).read ) ilist.erase(it++);
            else ++it;
        } 
    }
    for ( size_t i = 0; i < rids.size(); i++ ) {
        for ( it = dlist.begin(); it != dlist.end(); ) {
            if ( rids[i] == (*it).read ) dlist.erase(it++);
            else ++it;
        } 
    }

    __listToArray(ilist, INSERTION);
    __listToArray(dlist, DELETION);
}

void SpaPath::appendIndels( std::list<Mismatch> &mm, int type )
{
    if ( mm.size() == 0 ) return;

    std::list<Mismatch> mlist;
    std::list<Mismatch>::iterator it;

    if ( type == INSERTION ) mlist = std::list<Mismatch>( inss, inss+cins );
    else mlist = std::list<Mismatch>( dels, dels+cdel );
    
    for ( it = mm.begin(); it != mm.end(); ++it )
        mlist.push_back( *it );
    
    
    if ( type == INSERTION ) __listToArray(mlist, INSERTION);
    else __listToArray(mlist, DELETION);

}

void SpaPath::updateReadPlacement( int index,
                                  int start,
                                  AlignSummary &summary,
                                  bool verbose )
{
    if ( verbose ) {
        std::cout << "\t#ins:" << summary.ins.size() << "\t";
        std::cout << "\t#del:" << summary.del.size() << "\t";
        std::cout << "\t#mat:" << summary.match << "\t";
        std::cout << "\t#pos:" << summary.positive << "\t";
        std::cout << "\t#mms:" << summary.mismatch << "\n";
        std::cout << "\tseq1:" << summary.s1se.first << "\t" << summary.s1se.second << "\n";
        std::cout << "\tseq2:" << summary.s2se.first << "\t" << summary.s2se.second << "\n";
        std::cout << "\tlgap:" << summary.lgap.first << "\t" << summary.lgap.second << "\n";
        std::cout << "\tegap:" << summary.egap.first << "\t" << summary.egap.second << "\n";
        std::cout << "\t#out:" << summary.outer.first << "\t" << summary.outer.second << "\n";
        std::cout << "\t#in:" << summary.range.first << "\t" << summary.range.second << "\n";
        std::cout << "\tlength:" << summary.length << "\n";
        std::cout << "\t#score:" << summary.posrate << "\n";
    }
    //int s = start + summary.range.first;
    //inits[index] = 
    //__adjustReadStart( index, s, summary, verbose );
    __adjustReadStart( index, start, summary, verbose );
    dropMismatchRead(index);


    std::list<Mismatch>ilist, dlist;
    AlignPosList::iterator it;
    for ( it = summary.ilist.begin(); it != summary.ilist.end(); ++it ) {
        ilist.push_back( Mismatch(reads[index], it->seq_pos, start + it->ref_pos) );
    }
    for ( it = summary.dlist.begin(); it != summary.dlist.end(); ++it ) {
        dlist.push_back( Mismatch(reads[index], it->seq_pos, start + it->ref_pos) );
    }
    
    appendIndels( ilist, INSERTION );
    appendIndels( dlist, DELETION );
}

void SpaPath::dropReads( std::list<int> &index )
{
    if ( index.size() == 0 ) return;

    //std::tr1::unordered_map<int, bool> bad;
    for ( std::list<int>::iterator it = index.begin(); it != index.end(); ++it ) {
        //bad.insert( std::pair<int, bool>( *it, true ) );
        dropMismatchRead( *it );
    }


    //std::cout << "Old good reads:" << nread << "\n";

    std::map<int, bool> good;
    for ( size_t i = 0; i < nread; i++ ) 
        good.insert( std::pair<int, bool>(i, true) );

    for ( std::list<int>::iterator it = index.begin(); it != index.end(); ++it )
        good[*it] = false;

    ReadIdList ngood;
    for ( std::map<int, bool>::iterator it = good.begin(); it != good.end(); ++it ) {
        if ( it->second == true ) 
            ngood.push_back( reads[it->first] );
    }
    
    nread = ngood.size();
    ReadId *oreads = reads;
    reads = new ReadId[nread];
    int i = 0;
    for ( ReadIdList::iterator it = ngood.begin(); it != ngood.end(); ++it ) {
        reads[i] = *it; i++;
    }
    delete[] oreads;


    std::list<int> igood;
    for ( std::map<int, bool>::iterator it = good.begin(); it != good.end(); ++it ) {
        if ( it->second == true ) 
            igood.push_back( inits[it->first] );
    }
    
    int *oinits = inits;
    inits = new int[nread];
    i = 0;
    for ( std::list<int>::iterator it = igood.begin(); it != igood.end(); ++it ) {
        inits[i] = *it; i++;
    }
    delete[] oinits;

    



    //std::cout << "Old good reads:" << nread << "\n";

//     //nreads = SequenceArray(ngood.begin(), ngood.end());

    
//     ReadIdArray good = ReadIdArray( reads, reads+nread );
//     for ( int i = (int)good.size()-1; i >=0; i-- ) {
//         if ( bad.find(i) == bad.end() ) continue;

//         ReadIdArray temp;
//         int lsize, rsize;
//         rsize = good.size() - (i+1);  
//         lsize = i; 

//         if ( lsize > 0 && rsize > 0 ) {
//             temp.reserve(lsize+rsize );
//         } else if ( lsize > 0 ) {
//             temp.reserve(lsize );
//         } else if ( rsize > 0 ) {
//             temp.reserve(rsize );
//         }

//         if ( lsize > 0 ) {
//             temp.insert( temp.end(), good.begin(), good.begin() + lsize ); 
//         } else if ( rsize > 0 ) {
//             temp.insert( temp.end(), good.begin() + (i+1), good.end() ); 
//         }
//         good = temp;
//     }

//     nread = good.size();
//     std::cout << "New good reads:" << nread << "\n";
//     ReadId *old = reads;
//     reads = new ReadId[nread];
//     for ( unsigned i = 0; i < nread; i++ )
//         reads[i] = good[i];
    
//     delete[] old;
}

std::list<int> SpaPath::getPartialMatches(BitString *bstrs, double read_align_ratio, bool verbose)
{
//     if ( verbose ) {
//         std::cout << "Trimming partial aligned reads ...\n";
//         std::cout << consensus << "\n";
//         std::cout << "Consensus length:" << strlen(consensus) << "\n";
//     }

    std::list<int> partials;
    for ( size_t i = 0; i < nread; i++ ) {
        std::string query  = bstrs[reads[i]].toString();
        int seqlen = (int)query.size();

        int alnlen = seqlen;
        if ( inits[i] < 0 ) alnlen = seqlen + inits[i];
        else if ( inits[i]+seqlen >= (int)strlen(consensus) ) alnlen = strlen(consensus) - inits[i];

        if ( seqlen * read_align_ratio > alnlen ) {
            if ( verbose ) std::cout << "\tShort read:" << reads[i] << "\t" << query << "\tInit:" << inits[i] << "\tSeqLen:" << seqlen << "\tAlnLen:" << alnlen << "\n";
            partials.push_back(i); 
        }
    }

    return partials;
}

// void SpaPath::trimShortMatches(BitString *bstrs, double pcut, bool verbose)
// {
//     uList mis;
//     for ( unsigned i = 0; i < nread; i++ ) {
//         std::string query  = bstrs[reads[i]].toString();
//         int seqlen = (int)query.size();

//         int alnlen = seqlen;
//         if ( inits[i] < 0 ) alnlen = seqlen + inits[i];
//         else if ( inits[i]+seqlen >= strlen(consensus) ) alnlen = strlen(consensus) - inits[i];

//         if ( seqlen * pcut > alnlen ) {
//             if ( verbose ) std::cout << "\tshort:" << reads[i] << "\t" << query << "\t" << inits[i] << "\t" << alnlen << "\n";
//             mis.push_back(i); 
//         }
//     }
//     if ( verbose ) std::cout << "# Partial aligned reads:" << mis.size() << "\n";
//     __trim(mis, verbose);
// }

// Reset all reads to certain path 
void SpaPath::resetReads( PathId *used_reads,
                         PathId pid )
{
    for ( size_t i = 0; i < nread; i++ ) 
        used_reads[reads[i]] = pid;
}

void SpaPath::makeMismatchMap( std::tr1::unordered_map<ReadId, std::list<int> > & mmmap, int type )
{
    int num;
    Mismatch *mms;
    if ( type == DELETION ) {
        num = cdel; mms = dels;
    } else {
        num = cins; mms = inss;
    }

    for ( int i = 0; i < num; i++ ) {
        Mismatch mm = mms[i];
        if ( mmmap.find(mm.read) == mmmap.end()) 
            mmmap[mm.read] = std::list<int>();
        mmmap[mm.read].push_back(mm.spos);
    }
}

// void SpaPath::makeDeletionMap( std::tr1::unordered_map<ReadId, std::list<int> > & delmap )
// {
//     for ( unsigned i = 0; i < cdel; i++ ) {
//         Mismatch del = dels[i];
//         if ( delmap.find(del.read) == delmap.end()) 
//             delmap[del.read] = std::list<int>();
//         delmap[del.read].push_back(del.spos);
//     }

// }


// TO DO:
// Make full padded sueqneces.
// Currently, insertion AA is removed.
void SpaPath::setPaddedReads( std::vector<std::string> &nstrs, BitString *bstrs )
{
    std::tr1::unordered_map<ReadId, std::list<int> > delmap, insmap;
    //makeDeletionMap(delmap);
    makeMismatchMap(delmap, DELETION);
    makeMismatchMap(insmap, INSERTION);

    std::list<std::string> lreads;
    for ( unsigned i = 0; i < nread; i++ ) {
        std::string padded = std::string( strlen(consensus), '.' );
        std::string rstr = bstrs[reads[i]].toString();
        //size_t olen = rstr.length();
        std::list<int>::reverse_iterator it;
        for ( it = delmap[reads[i]].rbegin(); it != delmap[reads[i]].rend(); ++it ) {
            assert( *it >= 0 && *it < (int)rstr.size() );
            rstr.insert(*it, 1, '-');
        }
        for ( it = insmap[reads[i]].rbegin(); it != insmap[reads[i]].rend(); ++it ) {
            assert( *it >= 0 && *it < (int)rstr.size() ) ;
            rstr.erase(*it, 1);
        }
        
        int start = inits[i];
        if ( inits[i] < 0 ) {
//             if ( (int)rstr.size() > -1*inits[i] ) {
//                 std::cout << "Something wrong\n";
//                 std::cout << "Read:" << reads[i] << "\t" << rstr << "\n";
//                 std::cout << "Init:" << inits[i] << "\n";
//             }
            assert( (int)rstr.size() >= -1*inits[i] );
            rstr.erase(0, -1*inits[i]);
            start = 0;
        }

        //std::cout << reads[i] << "\tstart:" << start << "\tpadded:" << padded.size() << "\told-len:" << olen << "\trstr-len:" << rstr.length() << "\t" << rstr << "\n";
        if ( padded.size() < start+rstr.length() ) {
            assert ( rstr.size()-(start+rstr.length()-padded.size()) > 0 );
            rstr = rstr.substr(0, rstr.size()-(start+rstr.length()-padded.size()));
        }
        assert( padded.size() >= start + rstr.length() );
        padded.replace( start, rstr.length(), rstr);
        lreads.push_back(padded);
    }
    nstrs = std::vector<std::string>(lreads.begin(), lreads.end());
}

/*
 * Caution: insertion AA in tempoararily removed from a read sequence.
 */
void SpaPath::printAlignment(std::ostream &out, int csize, BitString *bstrs)
{
    std::string constr = std::string(consensus);
    if ( constr.size() == 0 ) return;

    std::vector<std::string> padded;
    setPaddedReads(padded, bstrs);

    std::multimap<int, unsigned> imap;
    for ( unsigned i = 0; i < nread; i++ )
        imap.insert( std::pair<int, unsigned>( inits[i], i ) );

    int from = 0;
    char buf[25];
    while ( from < (int)constr.size()-1 ) {
        if ( from+csize > (int)constr.length() )
            csize = constr.length()-from;

        assert(from+csize <= (int)constr.length());
        std::string rstr = constr.substr(from, csize);

        sprintf(buf, "\n%10d", from);
        out << buf << "\t";
        for ( size_t i = 0; i < rstr.size(); i++ ) {
            if ( (i+1)%5 == 0 && (i+1)%10 != 0 ) out << ".";
            else if ( (i+1)%10 == 0 ) out << ":";
            else out << " ";
        }
        out << "\n";

        sprintf(buf, "%10s", "consensus");
        out << buf << "\t" << rstr << "\n";
        for ( std::multimap<int, unsigned>::iterator it = imap.begin(); it != imap.end(); ++it ) {
            std::string read = padded[it->second];
            assert(from+csize <= (int)read.size());
            std::string sstr = read.substr(from, csize);
            if ( sstr != std::string(sstr.size(), '.') ) {
                sprintf(buf, "%10d", reads[it->second]);
                out << buf << "\t" << sstr << "\n"; 
            }
        }
        from += csize;
    }
}

void SpaPath::resetProfile()
{
    vprof->clear();
}

void SpaPath::setProfile(Profile &profile)
{
    vprof->build(profile);
}

void SpaPath::resetMismatches()
{
    if ( cins > 0  ) delete[] inss;
    if ( cdel > 0  ) delete[] dels;
    cins = cdel = 0;
}

void SpaPath::setMismatches( std::list<Mismatch> &mlist, int type )
{
    __listToArray( mlist, type );
}
