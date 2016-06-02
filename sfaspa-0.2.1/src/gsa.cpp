#include "gsa.h"

GSA::GSA()
{
    init(NULL, 0);
}

GSA::GSA( char **s, int n )
{
    init(s, n);
    GSA::buildSFA();
}

GSA::GSA( char **s, int n, bool gen_lcp )
{
    init(s, n);
    GSA::buildSFA();
    GSA::buildLCPs();
}

GSA::GSA( char **s, 
          int n, 
          bool gen_lcp, 
          const char *gsa_file,
          const char *lcp_file,
          const char *mcp_file)
{
    init(s,n);
    GSA::buildSFA( gsa_file );
    GSA::buildLCPs( lcp_file, mcp_file );
}

GSA::~GSA()
{
    clear();
}

void GSA::init(char **s, int n)
{
    Ids    = NULL;
    Pos    = NULL;
    seqs   = s;
    nreads = n;
}

void GSA::clear()
{
    if ( Ids != NULL ) delete[] Ids;
    if ( Pos != NULL ) delete[] Pos;
    Ids = NULL;
    Pos = NULL;
}

void GSA::buildSFA()
{
    double t0 = mytime();
    SFA::buildSFA();
    if ( verbose ) std::cout << "SFA created:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    //__convertWithBinarySearch();
    __convertWithArrays();
    if ( verbose ) std::cout << "GSA created:" << mytime()-t0 << " sec\n";
}

void GSA::buildSFA(const char *gsa_file)
{
    double t0 = mytime();
    SFA::buildSFA();
    if ( verbose ) std::cout << "SFA created:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    //__convertWithBinarySearch();
    __convertWithArrays();
    if ( verbose ) std::cout << "GSA created:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    GSA::writeSFA( gsa_file );
    GSA::clear(); 
    if ( verbose ) std::cout << "GSA written & purged:" << mytime()-t0 << " sec\n";
}

void GSA::buildLCPs()
{
    double t0 = mytime();
    SFA::buildLCP();
    if ( verbose ) std::cout << "LCP created:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    SFA::buildMLCP();
    if ( verbose ) std::cout << "lLCP/rLCP created:" << mytime()-t0 << " sec\n";
}

void GSA::buildLCPs( const char *lcp_file, const char *mcp_file )
{
    double t0 = mytime();
    SFA::buildLCP();
    if ( verbose ) std::cout << "LCP created:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    SFA::writeLCP( lcp_file );
    if ( verbose ) std::cout << "LCP written & purged:" << mytime()-t0 << " sec\n";
    t0 = mytime();
    SFA::purgeDoc();
    SFA::purgeSA();
    if ( verbose ) std::cout << "Concatenated string & Suffix array written/purged:" << mytime()-t0 << " sec\n";

    t0 = mytime();
    SFA::buildMLCP();
    if ( verbose ) std::cout << "lLCP/rLCP created:" << mytime()-t0 << " sec\n";
    t0 = mytime();
    SFA::writeMLCP( mcp_file );
    SFA::purgeMLCP();
    if ( verbose ) std::cout << "lLCP/rLCP written & purged:" << mytime()-t0 << " sec\n";
}

void GSA::__convertWithBinarySearch()
{    
    std::vector<SfaType> poss(nreads,0);
    size_t doc = 0;
    for ( size_t i = 0; i < size; i++ ) {
        if ( concat[i] == '$' ) {
            poss[doc] = i; doc++; 
        }
    }
 
    Ids = new SfaType[size];
    Pos = new LcpType[size];
    SfaType pos = 0;
    std::vector<SfaType>::iterator up;
    for ( size_t i = 0; i < size; i++ ) {
        pos = SA[i];
        up = std::lower_bound( poss.begin(), poss.end(), pos );
        doc = up-poss.begin();
        if ( doc > 0 ) 
            pos -= (poss[doc-1]+1);
        Ids[i] = doc;
        Pos[i] = pos;
    }
}

void GSA::__convertWithArrays()
{
    std::vector<SfaType> docs(size,0);
    std::vector<LcpType> poss(size,0);

    SfaType doc = 0;
    LcpType pos = 0;
    for ( size_t i = 0; i < size; i++ ) {
        docs[i] = doc;
        poss[i] = pos;
        pos++;
        if ( concat[i] == '$' ) {
            doc++; pos = 0;
        }
    }

    Ids = new SfaType[size];
    Pos = new LcpType[size];
    for ( size_t i = 0; i < size; i++ ) {
        SfaType doc = docs[SA[i]];
        LcpType pos = poss[SA[i]];
        Ids[i] = doc;
        Pos[i] = pos;
    }
}

void GSA::printSFA() 
{
    for ( size_t i = 0; i < size; i++ ) 
        std::cout << Ids[i] << ":" << (int)Pos[i] << " ";
    std::cout << "\n";
}

void GSA::printSuffix()
{
    // no reference sequence is given
    if ( seqs == NULL ) 
        return;

    for ( size_t i = 0; i < size; i++ ) {
        int doc = Ids[i];
        int pos = Pos[i];
        int rlen = strlen(seqs[doc]);
        assert(pos <= rlen);

        std::cout << i << "\t" << doc << "\t" << pos << "\t" << seqs[doc]+pos << "\n";
    }
}

void GSA::readSFA( const char *filename )
{
    std::fstream in;
    fio::openFile( in, filename, std::ios::in | std::ios::binary );

    in.read((char*)&size, sizeof(size_t));
    assert(size>0);

    Ids = new SfaType[size];
    Pos = new LcpType[size];
    for ( size_t i = 0; i < size; i++ ) {
        in.read((char*)&(Ids[i]), sizeof(SfaType));
        in.read((char*)&(Pos[i]), sizeof(LcpType));
    }
    in.close();
}

void GSA::load( const char *lcp_file,
                const char *mcp_file,
                const char *gsa_file )
{
    double t0 = mytime();
    GSA::readSFA( gsa_file );
    if ( verbose ) std::cout << "Suffix array loaded:" << mytime()-t0 << "\n";
    t0 = mytime();
    SFA::readLCP( lcp_file );
    if ( verbose ) std::cout << "LCP array loaded:" << mytime()-t0 << "\n";
    t0 = mytime();
    SFA::readMLCP( mcp_file );
    if ( verbose ) std::cout << "lLCP/rLCP array loaded:" << mytime()-t0 << "\n";
}

void GSA::load( const char *sfa_file, 
                const char *doc_file,
                const char *lcp_file,
                const char *mcp_file,
                const char *gsa_file )
{
    SFA::load( sfa_file, doc_file, lcp_file, mcp_file );
    GSA::readSFA( gsa_file );
}

void GSA::writeSFA( const char *gsa_file )
{
    std::fstream out;
    fio::openFile( out, gsa_file, std::ios::out | std::ios::binary );
    
    out.write((char*)&size, sizeof(size_t));
    for ( size_t i = 0; i < size; i++ ) {
        out.write((char*)&Ids[i], sizeof(SfaType));
        out.write((char*)&Pos[i], sizeof(LcpType));
    }
    out.close();
}

void GSA::dump( const char *sfa_file, 
                const char *doc_file,
                const char *lcp_file,
                const char *mcp_file,
                const char *gsa_file )
{
    SFA::dump( sfa_file, doc_file, lcp_file, mcp_file );
    GSA::writeSFA( gsa_file );
}

char* GSA::getSuffix( size_t p ) const
{
    assert( p < size );
    if ( seqs == NULL ) return NULL;
    SfaType rid = Ids[p];
    LcpType pos = Pos[p];
    return seqs[rid]+pos;
}

char* GSA::getEntireSuffix( size_t p )
{
    assert( p < size );
    if ( seqs == NULL ) return NULL;
    return SFA::getSuffix(p);
}


GsaType GSA::getAt( size_t p ) const
{
    if ( p >= size ) {
        std::cout << "size error:" << size << "\trequest:" << p << "\n";
        exit(-1);
    }
    assert( p < size );
    return GsaType( Ids[p], Pos[p] );
}

BoundType GSA::search( const SfaChar *srch, int len )
{
    assert(seqs!=NULL);
    return ( LCP != NULL && mLCP != NULL ) ?
        this->searchWithLCPs( srch, len ) :
        this->searchOnGSA( srch, len, 1, size ) ;
}


BoundType GSA::searchWithLCPs( const SfaChar *srch, int len )
{
    assert(seqs!=NULL);
    SfaType left  = this->getLeftBoundWithLCPs( srch, len );
    SfaType right = this->getRightBoundWithLCPs( srch, len );

    return BoundType(left-1, right-1);
}

BoundType GSA::searchOnGSA( const SfaChar *srch, int len, SfaType min, SfaType max )
{
    assert(seqs!=NULL);
    SfaType left  = this->getLeftBoundOnGSA( srch, len, min, max );
    SfaType right = this->getRightBoundOnGSA( srch, len, min, max );

    return BoundType(left-1, right-1);
}

char* GSA::getSeq( SfaType p )
{
    assert((size_t)p < size);                                                  

    int doc = Ids[p];
    int pos = Pos[p];
    int rlen = strlen(seqs[doc]);
    assert(pos <= rlen);
    
    return seqs[doc]+pos;
}

SfaType GSA::getLeftBoundWithLCPs(const SfaChar* pat, int len)
{
    if ( verbose ) std::cout << "Left bound:\n";
	SfaType lLcp = lcp( (SfaChar*)getSeq(0), pat );
	SfaType rLcp = lcp( (SfaChar*)getSeq(size-1), pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= getSeq(0)[lLcp]) return 1;
    if (rLcp <  len && pat[rLcp] >= getSeq(size-1)[rLcp]) return size+1;

	size_t left  = 1;
	size_t right = size;
    SfaType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        size_t mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);

        if ( lLcp >= rLcp ) {
            lMcp = getMcp(left, mid);
            if ( lMcp >= lLcp ) 
                mLcp = lLcp + lcp( (SfaChar*)getSeq(mid-1)+lLcp, pat+lLcp );
            else {
                mLcp = lMcp;
            }
        } else {
            rMcp = getMcp(mid,right);
            if ( rMcp >= rLcp ) 
                mLcp = rLcp + lcp( (SfaChar*)getSeq(mid-1)+rLcp, pat+rLcp );
            else {
                mLcp = rMcp;
            }
        }

		if (mLcp >= len || pat[mLcp] <= getSeq(mid-1)[mLcp]) { 
            right = mid; rLcp = mLcp; 
        }
		else { 
            left  = mid; lLcp = mLcp; 
        }
        if (verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}


SfaType GSA::getRightBoundWithLCPs(const SfaChar* pat, int len)
{
    if ( verbose ) std::cout << "Right bound:\n";
	SfaType lLcp = lcp((SfaChar*)getSeq(0), pat);
	SfaType rLcp = lcp((SfaChar*)getSeq(size-1), pat);
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= getSeq(0)[lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= getSeq(size-1)[rLcp]) return size;

	size_t left = 1;
	size_t right = size;
    SfaType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        size_t mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);
        if ( rLcp >= lLcp ) {
            rMcp = getMcp(mid,right);
            if ( rMcp >= rLcp ) 
                mLcp = rLcp + lcp( (SfaChar*)getSeq(mid-1)+rLcp, pat+rLcp );
            else {
                mLcp = rMcp;
            }
        }
        else {
            lMcp = getMcp(left, mid);
            if ( lMcp >= lLcp ) 
                mLcp = lLcp + lcp( (SfaChar*)getSeq(mid-1)+lLcp, pat+lLcp );
            else {
                mLcp = lMcp;
            }
        } 
		if (mLcp >= len || pat[mLcp] >= getSeq(mid-1)[mLcp]) { 
            left  = mid; lLcp = mLcp; 
        }
		else { 
            right = mid; rLcp = mLcp; 
        }
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}


SfaType GSA::getLeftBoundOnGSA(const SfaChar* pat, int len, SfaType l, SfaType r)
{
    if ( l > r ) return l;

    if ( verbose ) std::cout << "Left bound:\n";
	int lLcp = lcp( (SfaChar*)getSeq(0), pat );
	int rLcp = lcp( (SfaChar*)getSeq(size-1), pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= getSeq(0)[lLcp]) return 1;
    if (rLcp <  len && pat[rLcp] >= getSeq(size-1)[rLcp]) return size+1;

	// SfaType left = 1;
	// SfaType right = size;
	size_t left = l;
	size_t right = r;
	while (right-left > 1)
	{
		//size_t mid = left + floor((right-left)/2.0);
        size_t mid = (right+left)/2.0;
		int mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp( (SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
		if (mLcp >= len || pat[mLcp] <= getSeq(mid-1)[mLcp]) { 
            right = mid; rLcp = mLcp; 
        }
		else { 
            left  = mid; lLcp = mLcp; 
        }
        if (verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

// SfaType GSA::getLeftBoundOnGSA(const SfaChar* pat, int len, SfaType lmin, SfaType rmax)
// {
//     if ( rmax < lmin ) return lmin;
//     if ( verbose ) std::cout << "Left bound:\n";
// 	int lLcp = lcp( (SfaChar*)getSeq(lmin), pat );
// 	int rLcp = lcp( (SfaChar*)getSeq(rmax), pat );
//     if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

// 	if (lLcp >= len || pat[lLcp] <= getSeq(lmin)[lLcp]) return lmin+1;
//     if (rLcp <  len && pat[rLcp] >= getSeq(rmax)[rLcp]) return rmax+2;

// 	SfaType left = lmin+1;
// 	SfaType right = rmax+1;
// 	while (right-left > 1)
// 	{
// 		SfaType mid = left + floor((right-left)/2.0);
// 		int mLcp = lLcp <= rLcp? lLcp : rLcp;
// 		mLcp += lcp( (SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
//         if ( verbose ) 
//             printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
// 		if (mLcp >= len || pat[mLcp] <= getSeq(mid-1)[mLcp]) { 
//             right = mid; rLcp = mLcp; 
//         }
// 		else { 
//             left  = mid; lLcp = mLcp; 
//         }
//         if (verbose ) 
//             printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
// 	}

//     if ( verbose ) std::cout << "left bound:" << right << "\n";
// 	return right;
// }

SfaType GSA::getRightBoundOnGSA(const SfaChar* pat, int len, SfaType l, SfaType r) 
{
    if ( l > r ) return r;

    if ( verbose ) std::cout << "Right bound:\n";
	int lLcp = lcp((SfaChar*)getSeq(0), pat);
	int rLcp = lcp((SfaChar*)getSeq(size-1), pat);
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= getSeq(0)[lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= getSeq(size-1)[rLcp]) return size;

	// SfaType left = 1;
	// SfaType right = size;
	size_t left = l;
	size_t right = r;
	while (right-left > 1)
	{
		//size_t mid = left + floor((right-left)/2.0);
        size_t mid = (right+left)/2.0;
		int mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp((SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
		if (mLcp >= len || pat[mLcp] >= getSeq(mid-1)[mLcp]) { 
            left  = mid; lLcp = mLcp; 
        }
		else { 
            right = mid; rLcp = mLcp; 
        }
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

// SfaType GSA::getRightBoundOnGSA(const SfaChar* pat, int len, SfaType lmin, SfaType rmax) 
// {
//     if ( lmin > rmax ) return rmax;
//     if ( verbose ) std::cout << "Right bound:\n";
// 	int lLcp = lcp((SfaChar*)getSeq(lmin), pat);
// 	int rLcp = lcp((SfaChar*)getSeq(rmax), pat);
//     if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

// 	if (lLcp <  len && pat[lLcp] <= getSeq(lmin)[lLcp]) return lmin;
// 	if (rLcp >= len || pat[rLcp] >= getSeq(rmax)[rLcp]) return rmax+1;

// 	SfaType left = lmin+1;
// 	SfaType right = rmax+1;
// 	while (right-left > 1)
// 	{
// 		SfaType mid = left + floor((right-left)/2.0);
// 		int mLcp = lLcp <= rLcp? lLcp : rLcp;
// 		mLcp += lcp((SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
//         if ( verbose ) 
//             printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
// 		if (mLcp >= len || pat[mLcp] >= getSeq(mid-1)[mLcp]) { 
//             left  = mid; lLcp = mLcp; 
//         }
// 		else { 
//             right = mid; rLcp = mLcp; 
//         }
//         if ( verbose ) 
//             printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

// 	}

//     if ( verbose ) std::cout << "right bound:" << left << "\n";
// 	return left;
// }



BoundType GSA::refine( const SfaChar *srch, int len, SfaType lmin, SfaType rmax )
{
    SfaType left  = this->refineLeftBound( srch, len, lmin+1, rmax+1 );
    SfaType right = this->refineRightBound( srch, len, lmin+1, rmax+1 );

    return BoundType(left-1, right-1);
    
    //return this->searchOnGSA(srch, len, lmin+1, rmax+1);
}


SfaType GSA::refineLeftBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax)
{
    if ( lmin < 0 ) return lmin;
    if ( rmax < lmin ) return lmin;
    if ( verbose ) std::cout << "Left bound:" << pat << "\n";

    char *lstr = getSeq(lmin-1);
    char lch = ( (int)strlen(lstr) >= len ) ? lstr[len-1] : '\0';
	// char lch = getSeq(lmin-1)[len-1];
    char *rstr = getSeq(rmax-1);
    char rch = ( (int)strlen(rstr) >= len ) ? rstr[len-1] : '\0';

	// char lch = getSeq(lmin-1)[len-1];
	// char rch = getSeq(rmax-1)[len-1];
    char pch = pat[len-1];
    if ( verbose ) printf("lmin:%d\trmax:%d\tlch:%c\trch:%c\tpat:%c\n", lmin, rmax, lch, rch, pch);

    if ( pch >  rch ) return rmax+1;
    if ( pch <  lch ) return lmin+1;
    if ( pch == lch ) return lmin;
    
	SfaType left = lmin;
	SfaType right = rmax;
	while (right-left > 1)
	{
		SfaType mid = left + floor((right-left)/2.0);
		//char mch = lch <= rch? lch : rch;
        char *mstr = getSeq(mid-1);
        char mch = ( (int)strlen(mstr) >= len ) ? mstr[len-1] : '\0';
		//char mch = getSeq(mid-1)[len-1];
        if ( verbose ) 
            printf("left:%d right:%d mid:%d lch:%c rch:%c\n", left, right, mid, lch, rch);
		if ( pch <= mch ) {
            right = mid; rch = mch;
        }
		else { 
            left  = mid; lch = mch;
        }
        if (verbose ) 
            printf("left:%d right:%d mid:%d lch:%c rch:%c mch:%c\n", left, right, mid, lch, rch, mch );
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

SfaType GSA::refineRightBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax)
{
    if ( rmax < lmin ) return rmax;
    if ( verbose ) std::cout << "Right bound:" << pat << "\n";

    char *lstr = getSeq(lmin-1);
    char lch = ( (int)strlen(lstr) >= len ) ? lstr[len-1] : '\0';
	// char lch = getSeq(lmin-1)[len-1];
    char *rstr = getSeq(rmax-1);
    char rch = ( (int)strlen(rstr) >= len ) ? rstr[len-1] : '\0';

	//char rch = getSeq(rmax-1)[len-1];
    char pch = pat[len-1];
    if ( verbose ) printf("lmin:%d\trmax:%d\tlch:%c\trch:%c\tpat:%c\n", lmin, rmax, lch, rch, pch);

    if ( pch >  rch ) return rmax;
    if ( pch <  lch ) return lmin;
    if ( pch == rch ) return rmax;

	SfaType left = lmin;
	SfaType right = rmax;
	while (right-left > 1)
	{
		SfaType mid = left + floor((right-left)/2.0);
		//char mch = lch <= rch? lch : rch;
		//char mch = getSeq(mid-1)[len-1];
        char *mstr = getSeq(mid-1);
        char mch = ( (int)strlen(mstr) >= len ) ? mstr[len-1] : '\0';

        if ( verbose ) 
            printf("left:%d right:%d mid:%d lch:%c rch:%c\n", left, right, mid, lch, rch);
		if (  pch >= mch ) {
            left  = mid; lch = mch; 
        }
		else { 
            right = mid; rch = mch; 
        }
        if (verbose ) 
            printf("left:%d right:%d mid:%d lch:%c rch:%c mch:%c\n", left, right, mid, lch, rch, mch );
	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

BoundType GSA::getEndBound( const SfaChar *srch, int len, SfaType lmin, SfaType rmax )
{
    SfaType left  = this->getLeftEndBound( srch, len, lmin+1, rmax+1 );
    SfaType right = this->getRightEndBound( srch, len, lmin+1, rmax+1 );

    return BoundType(left-1, right-1);
}


SfaType GSA::getLeftEndBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax)
{
    //if ( rmax < lmin ) return lmin+1;
    //if ( verbose ) std::cout << "Left bound:" << pat << "\n";
    if ( lmin < 0 ) return lmin;
    if ( rmax < lmin ) return lmin;

    assert( lmin >= 1 );
    assert( (size_t)rmax <= size);

    const char *lstr = getSeq(lmin-1);
    const char *rstr = getSeq(rmax-1);
    
    if ( strlen(rstr)+1 < (size_t)len ) return rmax+1;

	// char lch = getSeq(lmin-1)[len];
	// char rch = getSeq(rmax-1)[len];
    char lch = ( strlen(lstr)+1 >= (size_t)len ) ? lstr[len] : 127;
    char rch = rstr[len];
    char pch = pat[len];

    if ( lch == '\0' ) return lmin;
    if ( lch > '\0' ) return 0;
    
    if ( verbose ) printf("lmin:%d\trmax:%d\tlch:%c\trch:%c\tpat:%c\n", lmin, rmax, lch, rch, pch);

	size_t left = lmin;
	size_t right = rmax;
	while (right-left > 1)
	{
		size_t mid = left + floor((right-left)/2.0);
        assert( mid>=0 && mid<size);
		//char mch = lch <= rch? lch : rch;
        const char *mstr = getSeq(mid-1);
        size_t mlen = strlen(mstr);
		//char mch = getSeq(mid-1)[len];
        char mch = ( mlen+1 >= (size_t)len ) ? mstr[len] : 127;
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lch:%c rch:%c\n", left, right, mid, lch, rch);
		if ( mlen+1 >= (size_t)len && '\0' <= mch ) {
            right = mid; rch = mch;
        }
		else { 
            left  = mid; lch = mch;
        }
        if (verbose ) 
            printf("left:%zu right:%zu mid:%zu lch:%c rch:%c mch:%c\n", left, right, mid, lch, rch, mch );
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

SfaType GSA::getRightEndBound(const SfaChar* pat, int len, SfaType lmin, SfaType rmax)
{
    //if ( rmax < lmin ) return rmax+1;
    //if ( verbose ) std::cout << "Right bound:" << pat << "\n";

    if ( rmax < lmin ) return rmax;
    
    assert( lmin >= 1 );
    assert( (size_t)rmax <= size );

    const char *lstr = getSeq(lmin-1);
    const char *rstr = getSeq(rmax-1);
    
    if ( strlen(rstr)+1 < (size_t)len ) return rmax;
    
	// char lch = getSeq(lmin-1)[len];
	// char rch = getSeq(rmax-1)[len];
    char lch = ( strlen(lstr)+1 >= (size_t)len ) ? lstr[len] : 127;
    char rch = rstr[len];
    char pch = pat[len];

    if ( rch == '\0' ) return rmax;
    if ( lch > '\0' ) return -1;
    
    if ( verbose ) printf("lmin:%d\trmax:%d\tlch:%c\trch:%c\tpat:%c\n", lmin, rmax, lch, rch, pch);

	//if ( pat[len] <= lch ) return lmin;
    //if ( pat[len] >= rch ) return rmax+1;

	size_t left = lmin;
	size_t right = rmax;
	while (right-left > 1)
	{
		size_t mid = left + floor((right-left)/2.0);
        assert( mid>=0 && mid<size);
		//char mch = lch <= rch? lch : rch;
        const char *mstr = getSeq(mid-1);
        size_t mlen = strlen(mstr);
		//char mch = getSeq(mid-1)[len];
        char mch = ( mlen+1 >= (size_t)len ) ? mstr[len] : 127;
		// char mch = getSeq(mid-1)[len];
        if ( verbose ) 
            printf("left:%zu right:%zu mid:%zu lch:%c rch:%c\n", left, right, mid, lch, rch);
		//if (  '\0' >= mch ) {
		if ( mlen+1 >= (size_t)len && '\0' >= mch ) {
            left  = mid; lch = mch; 
        }
		else { 
            right = mid; rch = mch; 
        }
        if (verbose ) 
            printf("left:%zu right:%zu mid:%zu lch:%c rch:%c mch:%c\n", left, right, mid, lch, rch, mch );
	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

