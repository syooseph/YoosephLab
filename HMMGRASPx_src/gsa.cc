#include "gsa.h"

GSA::GSA()
{
    init(NULL, 0);
}

GSA::GSA( char **s, RidType n )
{
    init(s, n);
    GSA::buildSFA();
}

GSA::GSA( char **s, RidType n, bool gen_lcp )
{
    init(s, n);
    GSA::buildSFA();
    GSA::buildLCPs();
}

GSA::GSA( char **s, 
          RidType n, 
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

void GSA::init(char **s, RidType n)
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
    std::vector<IdxType> poss(nreads,0);
    IdxType doc = 0;
    for ( IdxType i = 0; i < size; i++ ) {
        if ( concat[i] == '$' ) {
            poss[doc] = i; doc++; 
        }
    }
 
    Ids = new RidType[size];
    Pos = new LcpType[size];
    RidType pos = 0;
    std::vector<IdxType>::iterator up;
    for ( IdxType i = 0; i < size; i++ ) {
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
    std::vector<RidType> docs(size,0);
    std::vector<LcpType> poss(size,0);

    RidType doc = 0;
    LcpType pos = 0;
    for ( IdxType i = 0; i < size; i++ ) {
        docs[i] = doc;
        poss[i] = pos;
        pos++;
        if ( concat[i] == '$' ) {
            doc++; pos = 0;
        }
    }

    Ids = new RidType[size];
    Pos = new LcpType[size];
    for ( IdxType i = 0; i < size; i++ ) {
        RidType doc = docs[SA[i]];
        LcpType pos = poss[SA[i]];
        Ids[i] = doc;
        Pos[i] = pos;
    }
}

void GSA::printSFA() 
{
    for ( IdxType i = 0; i < size; i++ ) 
        std::cout << Ids[i] << ":" << (int)Pos[i] << " ";
    std::cout << "\n";
}

void GSA::printSuffix()
{
    // no reference sequence is given
    if ( seqs == NULL ) 
        return;

    for ( IdxType i = 0; i < size; i++ ) {
        RidType doc = Ids[i];
        LcpType pos = Pos[i];
        LcpType rlen = (LcpType) strlen(seqs[doc]);
        assert(pos <= rlen);

        std::cout << i << "\t" << doc << "\t" << pos << "\t" << seqs[doc]+pos << "\n";
    }
}

void GSA::readSFA( const char *filename )
{
    std::fstream in;
    fio::openFile( in, filename, std::ios::in | std::ios::binary );

    in.read((char*)&size, sizeof(IdxType));
    assert(size>0);

    Ids = new RidType[size];
    Pos = new LcpType[size];
    for ( IdxType i = 0; i < size; i++ ) {
        in.read((char*)&(Ids[i]), sizeof(RidType));
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
    
    out.write((char*)&size, sizeof(IdxType));
    for ( IdxType i = 0; i < size; i++ ) {
        out.write((char*)&Ids[i], sizeof(RidType));
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

char* GSA::getSuffix( IdxType p )
{
    assert( p <= size );
    if ( seqs == NULL ) return NULL;
    RidType rid = Ids[p];
    LcpType pos = Pos[p];
    return seqs[rid]+pos;
}

char* GSA::getEntireSuffix( IdxType p )
{
    assert( p <= size );
    if ( seqs == NULL ) return NULL;
    return SFA::getSuffix(p);
}


GsaType GSA::getAt( IdxType p )
{
    assert( p <= size );
    return GsaType( Ids[p], Pos[p] );
}

BoundType GSA::search( const SfaChar *srch, LcpType len )
{
    assert(seqs!=NULL);
    return ( LCP != NULL && mLCP != NULL ) ?
        this->searchWithLCPs( srch, len ) :
        this->searchOnGSA( srch, len ) ;
}

BoundType GSA::searchWithLCPs( const SfaChar *srch, LcpType len )
{
    assert(seqs!=NULL);
    IdxType left  = this->getLeftBoundWithLCPs( srch, len );
    IdxType right = this->getRightBoundWithLCPs( srch, len );

    return BoundType(left-1, right-1);
}

BoundType GSA::searchWithLCPs_bounded(
    const SfaChar *srch, LcpType len, 
    IdxType left_bound, IdxType right_bound
) {
  assert(seqs!=NULL);
  IdxType left  = this->getLeftBoundWithLCPs_bounded(srch, len, left_bound, right_bound);
  IdxType right = this->getRightBoundWithLCPs_bounded(srch, len, left_bound, right_bound);
  return BoundType(left-1, right-1);
}

BoundType GSA::searchOnGSA( const SfaChar *srch, LcpType len )
{
    assert(seqs!=NULL);
    IdxType left  = this->getLeftBoundOnGSA( srch, len );
    IdxType right = this->getRightBoundOnGSA( srch, len );

    return BoundType(left-1, right-1);
}

char* GSA::getSeq(IdxType p )
{
    assert(p < size);                                                  

    RidType doc = Ids[p];
    LcpType pos = Pos[p];
    LcpType rlen = (LcpType) strlen(seqs[doc]);
    assert(pos <= rlen);
    
    return seqs[doc]+pos;
}

IdxType GSA::getLeftBoundWithLCPs(const SfaChar* pat, LcpType len)
{
    if ( verbose ) std::cout << "Left bound:\n";
	  LcpType lLcp = lcp( (SfaChar*)getSeq(0), pat );
	  LcpType rLcp = lcp( (SfaChar*)getSeq(size-1), pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= getSeq(0)[lLcp]) return 1;
    if (rLcp <  len && pat[rLcp] >= getSeq(size-1)[rLcp]) return size+1;

	IdxType left  = 1;
	IdxType right = size;
  LcpType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        IdxType mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);

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
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

IdxType GSA::getLeftBoundWithLCPs_bounded(
    const SfaChar* pat, LcpType len, 
	  IdxType left_bound, IdxType right_bound
) {
    if ( verbose ) std::cout << "Left bound:\n";
	LcpType lLcp = lcp( (SfaChar*)getSeq(0), pat );
	LcpType rLcp = lcp( (SfaChar*)getSeq(size-1), pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= getSeq(0)[lLcp]) return 1;
    if (rLcp <  len && pat[rLcp] >= getSeq(size-1)[rLcp]) return size+1;

	IdxType left  = left_bound + 1;
	IdxType right = right_bound + size;
  LcpType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        IdxType mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);

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
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

IdxType GSA::getRightBoundWithLCPs(const SfaChar* pat, LcpType len) 
{
    if ( verbose ) std::cout << "Right bound:\n";
	LcpType lLcp = lcp((SfaChar*)getSeq(0), pat);
	LcpType rLcp = lcp((SfaChar*)getSeq(size-1), pat);
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= getSeq(0)[lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= getSeq(size-1)[rLcp]) return size;

	size_t left = 1;
	size_t right = size;
    LcpType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        IdxType mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);
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
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

IdxType GSA::getRightBoundWithLCPs_bounded(
    const SfaChar* pat, LcpType len, 
	  IdxType left_bound, IdxType right_bound
) {
    if ( verbose ) std::cout << "Right bound:\n";
	LcpType lLcp = lcp((SfaChar*)getSeq(0), pat);
	LcpType rLcp = lcp((SfaChar*)getSeq(size-1), pat);
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= getSeq(0)[lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= getSeq(size-1)[rLcp]) return size;

	IdxType left = left_bound + 1;
	IdxType right = right_bound + 1;
    LcpType mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        IdxType mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);
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
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

IdxType GSA::getLeftBoundOnGSA(const SfaChar* pat, LcpType len)
{
    if ( verbose ) std::cout << "Left bound:\n";
	LcpType lLcp = lcp( (SfaChar*)getSeq(0), pat );
	LcpType rLcp = lcp( (SfaChar*)getSeq(size-1), pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= getSeq(0)[lLcp]) return 1;
    if (rLcp <  len && pat[rLcp] >= getSeq(size-1)[rLcp]) return size+1;

	IdxType left = 1;
	IdxType right = size;
	while (right-left > 1)
	{
		IdxType mid = left + floor((right-left)/2.0);
		LcpType mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp( (SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
        if ( verbose ) 
            printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
		if (mLcp >= len || pat[mLcp] <= getSeq(mid-1)[mLcp]) { 
            right = mid; rLcp = mLcp; 
        }
		else { 
            left  = mid; lLcp = mLcp; 
        }
        if (verbose ) 
            printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

IdxType GSA::getRightBoundOnGSA(const SfaChar* pat, LcpType len) 
{
    if ( verbose ) std::cout << "Right bound:\n";
	LcpType lLcp = lcp((SfaChar*)getSeq(0), pat);
	LcpType rLcp = lcp((SfaChar*)getSeq(size-1), pat);
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= getSeq(0)[lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= getSeq(size-1)[rLcp]) return size;

	IdxType left = 1;
	IdxType right = size;
	while (right-left > 1)
	{
		IdxType mid = left + floor((right-left)/2.0);
		LcpType mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp((SfaChar*)getSeq(mid-1)+mLcp, pat+mLcp);
        if ( verbose ) 
            printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d\n", left, right, mid, lLcp, rLcp);
		if (mLcp >= len || pat[mLcp] >= getSeq(mid-1)[mLcp]) { 
            left  = mid; lLcp = mLcp; 
        }
		else { 
            right = mid; rLcp = mLcp; 
        }
        if ( verbose ) 
            printf("left:%d right:%d mid:%d lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

//  ############################################################################
//  extended by Cuncong Zhong

LcpType  GSA::getSuffixLength( IdxType p )
{
    assert( p <= size );
    if ( seqs == NULL ) return 0;
    RidType rid = Ids[p];//GSFA[p].doc;
    LcpType pos = Pos[p];//GSFA[p].pos;
    return strlen(seqs[rid]+pos);
}

char* GSA::getSuffix_explicit(RidType rid, LcpType pos) {
  assert(rid >= 0 && rid < nreads);
  assert(pos >= 0);
  if ( seqs == NULL ) return NULL;
  return seqs[rid]+pos;
}

RidType GSA::getId(IdxType p)  {
  return Ids[p];
}

LcpType GSA::getPos(IdxType p) {
  return Pos[p];
}

LcpType GSA::getFullSequenceLength(IdxType p)  {
  assert(p >= 0 && p < size);
  return (LcpType) strlen(seqs[Ids[p]]);
}

char* GSA::getSequence_explicit(RidType rid) {
  assert(rid >= 0 && rid < nreads);
  if ( seqs == NULL ) return NULL;
  return seqs[rid];
}

LcpType GSA::getLcp(IdxType p)
{
    assert( p <= size );
    return LCP[p];
}

LcpType GSA::getSeqLength_RID(RidType rid) {
  return strlen(seqs[rid]);
}


