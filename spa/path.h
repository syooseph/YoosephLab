//==============================================================================
// Wed 2011-06-08 08:04:21 PM
//==============================================================================

#ifndef __SPAPATH_H__
#define __SPAPATH_H__

#include <iostream>
#include <list>
#include <set>
#include "kmer.h"
#include "InvertedIndex.h"
#include "set.h"
#include "biostr.h"
#include "bitstr.h"
#include "galignment.h"
#include "lalignment.h"
#include "param.h"
//#incldue "alpha.h"
#include "vprofile.h"

typedef int PathId;
const PathId NOT_PATH = -1;
const ReadId NOT_PAIR = 4294967295; // 2^32-1 

enum INDEL { INSERTION, DELETION };
enum DIRECTION { dLEFT, dRIGHT };

typedef std::list<unsigned> uList;
typedef std::pair<PathId, AlignSummary> PathAlignPair;
typedef std::list<PathAlignPair> PathAlignPairList; 
typedef std::tr1::unordered_map<PathId, PathAlignPairList> PathBins;

struct Mismatch {
	ReadId   read; // read id
	int spos;
	int rpos;
	
	Mismatch(){}

	Mismatch( ReadId rid, int sp, int rp ) 
	{
		read = rid; spos = sp; rpos = rp;
	}

	Mismatch( const Mismatch &source )
	{
		read = source.read; spos = source.spos; rpos = source.rpos;
	}

	Mismatch& operator= ( const Mismatch &source )
	{
		read = source.read; spos = source.spos; rpos = source.rpos;
		return *this;
	}

	void dump(std::ostream &out) 
	{
		out.write((char*)&read, sizeof(ReadId));
		out.write((char*)&spos, sizeof(int));
		out.write((char*)&rpos, sizeof(int));
	}

	void load(std::istream &in) 
	{
		in.read((char*)&read, sizeof(ReadId));
		in.read((char*)&spos, sizeof(int));
		in.read((char*)&rpos, sizeof(int));
	}
};

typedef std::pair<int, int> IntPair;


class SpaPath
{
 private:
	char*  consensus;
	unsigned   nread; // no. of reads
	unsigned   nkmer; // no. of kmers
	KmerId*    kmers; // kmers
	ReadId*    reads; // reads
	int*       inits; // alignment start pos of each read
	Mismatch*   inss;
	Mismatch*   dels;
	unsigned    cins; // no. of insertion reads
	unsigned    cdel; // no. of deleteion reads
	ProfileVector *vprof;

 public:
	SpaPath();
	SpaPath( KmerId *kids, ReadId *rids, unsigned nk, unsigned nr );
	SpaPath( const char *seq, KmerId *kids, ReadId *rids, unsigned nk, unsigned nr );
	SpaPath( const SpaPath &source );
	SpaPath& operator= ( const SpaPath &source );

	~SpaPath();
	void build(SpaPath&);
	void clear();
	void align(InvertedIndex &, BitString *bstrs, Param &param);
	void align(BitString *bstrs, Param &param);

	/* getters */
	char*       getConsensus()       { return consensus; }
	std::string getConsensusString() { return std::string(consensus); }
	unsigned    getReadCount()       { return nread; }
	unsigned    getKmerCount()       { return nkmer; }
	KmerId*     getKmers()           { return kmers; }
	ReadId*     getReads()           { return reads; }
	int*        getInits()           { return inits; }
	int         getInit(int i)       { return inits[i]; }
	Mismatch*   getInsertions()      { return inss; }
	Mismatch*   getDeletions()       { return dels; }
	unsigned    countInsertions()    { return cins; }
	unsigned    countDeletions()     { return cdel; }
	ProfileVector* getProfileVector() { return vprof; }

	/* setter */
	void setInit(unsigned i, int n) { inits[i] = n; }
	void setMismatches( std::list<Mismatch> &mlist, int type );

	/* other functions */
	void dump(std::ostream &);
	void load(std::istream &);

	void validate(BitString*, Param &param);
	
	void join( SpaPath *other, int start, int lgap, int tgap, AlignPosList &ilist, AlignPosList &dlist, int kmer_size,  InvertedIndex &iindex, BitString *bstrs, Param &param);
	void join( SpaPath *oaln, std::string &mid, int start, int lgap, int tgap, AlignPosList &ilist, AlignPosList &dlist, int kmer_size,  InvertedIndex &iindex, BitString *bstrs, Param &param );
	void addRead( ReadId rid, int spos );
	void addReads( ReadId *rids, int *sposs, size_t nrid );

	std::pair<int, PathId> getLeadingGapCount( PathAlignPairList &path_aligns );
	std::pair<int, PathId> getTrailingGapCount( PathAlignPairList &path_aligns );
	void mergeCluster( PathAlignPairList &path_aligns,
					   std::tr1::unordered_map<PathId, SpaPath*> &path2aln_map,
					   BitString *bstrs,
					   Param &param );
		

	void joinReads( ReadId *rids, int *poss, int rsize, BitString *bstrs, int k, int direction, bool verbose );
	void updateConsensus( std::string &str, BitString*, double, int k, bool );
	void updateConsensusSequence( const char * );

	void adjustInit(int i, int g);
	void updateKmers( KmerId *nkids, unsigned nsize );
	void adjustPositions(int);
	void adjustStartPositions(int);
	void trimIndels(size_t spos);

	void dropMismatchRead( int index );
	void dropMismatches( ReadIdArray & );
	void appendIndels( std::list<Mismatch> &mm, int type );

	void __adjustReadStart( int index, int init, AlignSummary &summary, bool verbose );
	void updateReadPlacement( int index, int start, AlignSummary &summary, bool verbose );
		
	void dropReads( std::list<int>& );

	//void trimShortMatches(BitString *bstrs, double pcut, bool verbose);
	std::list<int> getPartialMatches(BitString *bstrs, double , bool verbose);
	void resetReads( PathId *used_reads, PathId pid );
	
	void printAlignment(std::ostream &out, int csize, BitString *bstrs);

	void setProfile(Profile &profile);
	void resetProfile();

	void resetMismatches();

	void append(std::string &, size_t kmre_size);
	void prepend(std::string &that_seq, size_t kmer_size);
 private:
	void __copy( const SpaPath &source );
	void build(KmerId *kids, ReadId *rids, unsigned nk, unsigned nr);
	void resetIndels();
	void updateInsertionPos( AlignPosList &nlist );
	void makeGappyConsensus(AlignPosList &ninss, Param &param);
	bool __badAlignment( AlignSummary &summary, int qlen, Param &param );
	std::list<int> updateGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type, Param &param );
	void insertGapsToReads( AlignPosList &nlist, BitString *bstrs, int type );
	void insertGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type );
	void verifyGaps( AlignPosList &nlist, BitString *bstrs, IntPair range, int type );
	AlignSummary __realignRead( std::string &query, int start, int end, bool verbose );
	void __initFlags(bool**, InvertedIndex &);
	void __initFlags(bool **flags, BitString *bstrs, int kmer_size);
	void __destroy(bool**);
	void __trim(uList&, bool);
	uList __mismatch(bool**, unsigned, unsigned, unsigned, unsigned);
	unsigned __start(bool**, unsigned);
	unsigned __start(bool**, unsigned, unsigned);
	unsigned __end(bool**, unsigned);
	unsigned __end(bool**, unsigned, unsigned);
	unsigned __nstart( bool **flags, unsigned rid, unsigned s );
	std::string __alignedString(bool **, unsigned, unsigned, unsigned, unsigned);
	int __matchlength( std::string&, std::string& );
	/* unsigned __fitstart( std::string &, std::string &, bool ); */
	/* unsigned __alnStartPos(BitString*, unsigned, unsigned, unsigned, bool); */
	void __listToArray(std::list<Mismatch> &, int);
	IntPair __getRange( bool **flags, unsigned i, unsigned seq_len, int k );

	void mergeReads(SpaPath*);
	//void mergeIndel(AlignPosList&, AlignPosList&);
	void mergeIndels(SpaPath *oaln);
	void prepend(SpaPath*, int, int);
	void append(SpaPath*, int, int);

	KmerArray getKmerArray(std::string&, int);

	void __joinReadsLeft ( ReadId *rids, int *sposs, int rsize, BitString *bstrs, int k, bool verbose );
	void __joinReadsRight( ReadId *rids, int *rposs, int rsize, BitString *bstrs, int k, bool verbose );
	void __updateReadPos( std::string &, BitString*, double merge_score, bool verbose );

	void __updateIndel(unsigned i, int offset, AlignSummary &summary, std::list<Mismatch> &ilist, std::list<Mismatch> &dlist);
	bool __alignReadToPath(AlignSummary &summary, std::string &ref, unsigned i, bool **flags, BitString *bstrs, std::list<Mismatch> &ilist, std::list<Mismatch> &dlist, Param &param);
	bool __short( unsigned seq_len, unsigned aln_len, Param &param);
	bool __weak(IntPair &se, unsigned seq_len, unsigned aln_len, bool **flags, unsigned i, Param &param);
	AlignSummary __compareBase(std::string &query, std::string &sbjct);
	bool __doAlignment(AlignSummary &summary, std::string &query, std::string &sbjct, Param &param);

	void setPaddedReads( std::vector<std::string> &nstrs, BitString *bstrs );
	//void makeDeletionMap( std::tr1::unordered_map<ReadId, std::list<int> > & delmap );
	void makeMismatchMap( std::tr1::unordered_map<ReadId, std::list<int> > & mmmap, int type );

};


#endif

