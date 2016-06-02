//==============================================================================
// Thu 2011-06-16 10:21:48 AM
//==============================================================================

#ifndef __MSA_H__
#define __MSA_H__

#include "path.h"
#include "profile.h"
#include "param.h"

enum PLATFORM { ILLUMINA, SOLID, IONTORRENT, ROCHE, SANGER };

enum BUILDTASK { PILE = 1, FLAT = 2, TRIM = 4, CODON = 8 };

struct Task {
	bool pile, flat, trim, profile, codon;
	Task() { pile = flat = trim = profile = codon = false; }
};



typedef int SeqPos;
typedef int RefPos;
typedef unsigned Count;

typedef std::pair<bool, bool> FlagPair;
typedef std::map<RefPos, Count> PosCountsMap;
typedef std::tr1::unordered_map<ReadId, Count> ReadCountsMap;
typedef std::map<RefPos, ReadCountsMap> PosReadCountsMap;
typedef std::pair<SeqPos, RefPos> PosPair;
typedef std::list<PosPair> PosPairList;
typedef std::map<ReadId, PosPairList> ReadPosPairListMap;
typedef std::list<std::string> SequenceList;
typedef std::vector<std::string> SequenceArray;


struct InChar 
{
	int rpos;
	char ch;
	InChar( unsigned r, char c ) : rpos(r), ch(c) {}
};

/**
 * Multiple sequence alignment.
 */
class MSA
{
 private:
	std::string consensus;
	SequenceArray nreads;
	Profile profile;
	bool aligned;

	void init();
	void build(SpaPath &kaln, PathId pid, PathId *path_reads, BitString *bstrs, char *strands, ReadId *pairs, int mode, int &error, Param &param );
	Task parseMode( int mode );
	int minStartPosition( SpaPath &kaln );
	int maxEndPosition( SpaPath &kaln, BitString *bstrs, int k, bool verbose);
	void adjustStartPositions( SpaPath &kaln, int k, bool verbose );
	bool pileUpReads( SpaPath &kaln, BitString *bstrs, int k, bool verbose);
	FlagPair trimEnds( SpaPath &kaln, int min_depth, bool verbose );
	void adjustInsertions( SpaPath&, ReadPosPairListMap&, PosReadCountsMap&, PosCountsMap&, bool);
	void stretch( SpaPath &, BitString *bstrs, unsigned k, bool verbose);
	void stretchConsensus( SpaPath &kaln, BitString *bstrs, unsigned k, bool verbose);
	ReadPosPairListMap loadInsertions( SpaPath &, bool );
	PosReadCountsMap countReadPositions(ReadPosPairListMap &);
	PosCountsMap getMaxCounts(PosReadCountsMap &);
	void insertGapsToReads( SpaPath &, BitString*, bool);
	void insertGapsToRef(PosCountsMap &, SpaPath &, int, bool );
	void trimTail(bool);
	void fitLength(bool);
	void trimHead(SpaPath &, bool);
	void makeConsensus();
	void makeProfile();
	void updateAllPosFromConsensus( SpaPath &spath, bool verbose );
	void updateInit( SpaPath &spath, bool );
	void __initializeProfile();
	void __resetProfile();
	void __updateProfile();
	void __update(SpaPath &spath, int kmer_size, bool verbose);

	int __sumColumn( int col );
	int __countLeadingGaps(int);
	int __countTrailingGaps(int);

	void __trimReads(int gap, int direction);
	void __trimProfile(int gap, int direction);

	void __trimPoorRegion( std::vector<IntPair> &poors );
	std::vector<IntPair> __searchSimilarRegion( IntPair qr, bool verbose );
	std::vector<IntPair> __getSimilarRegions( std::string &query, int qs, int qe, int direction, bool verbose );
	bool __flatten( BitString *bstrs, char *strands, ReadId *pairs, SpaPath &spath, Param &param );
	bool __flattenSingleEnds( std::vector<IntPair> &poors, SpaPath &spath, BitString *bstrs, Param &param );
	bool __flattenPairedEnds( std::vector<IntPair> &poors, SpaPath &spath, BitString *bstrs, char *strands, ReadId *pairs, Param &param );
	std::vector<int> __getRandomSubset( std::vector<int> &index, size_t size, unsigned seed );
	std::vector<int> __getReadsInRegion( SpaPath &aln, BitString *bstrs, IntPair region, int offset );

	void __borrow( IntPair bad, std::vector<IntPair> &sim, SpaPath &kaln, BitString *bstrs, double merge_score, int offset, unsigned seed, bool verbose );
	void __placeToRegion( IntPair region, std::vector<int> &subset, SpaPath &kaln, BitString *bstrs, int offset, double merge_score, bool verbose );
	void __redistribute();

	bool __correctMisplacement(SpaPath &spath, BitString *bstrs, char *strands, ReadId *pairs, Param &param );
	int findPair( ReadId *reads, size_t nread, ReadId pread );

	void __resetRead(int index);	
	void __assignRead(int index, int start, std::string seq);

	void joinUsedReads( PathId pid, ReadId *reads, int nread, PathId *used_reads );
	void joinUsedReads( PathId pid, ReadIdArray &path_rids, PathId *used_reads );

	void dropReads( std::list<int>& );
	int __getMinAlignStart( SpaPath &kaln );

	void alignReads(  SpaPath &kaln, BitString *bstrs, int k, bool verbose);
	void updateConsensus(SpaPath &kaln, int k);

	bool trimAfterStopCodon( SpaPath &kaln, bool verbose );
	bool trimPartialMatches( SpaPath &kaln, BitString *bstrs, double merge_score, bool verbose );

 public:
	MSA();
	MSA(SpaPath &kaln, PathId pid, PathId *path_reads, BitString *bstrs, char *strands, ReadId *pairs, int mode, int &error, Param &param );
	~MSA();
	std::string getConsensus();
	Profile getProfile();
	void printProfile(std::ostream &);
	void printAlignment(std::ostream &, SpaPath &, int);
	std::vector<IntPair> getZeroCoverageRegions();
};

#endif


