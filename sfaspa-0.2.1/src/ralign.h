/**
 * \file      ralign.h
 * \brief     Aligned reads
 * \details   After graph path is extracted, reads are placed.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2013
 * \copyright J. Craig Venter Institute.
 */

#ifndef __READ_ALIGNER_H__
#define __READ_ALIGNER_H__

#include "core.h"
#include "log.h"
#include "rplace.h"
#include "rcmat.h"
#include "ReadStartCount.h"

typedef std::tr1::unordered_map<std::string, size_t> WordFreqMap;

/**
 * \brief Pair-wise alignment range between a path and a read
 */
struct AlignRange
{
	int path_beg; ///< alignment start position in a path
	int path_end; ///< alignment end position in a path
	int read_beg; ///< alignment start position in a read
	int read_end; ///< alignment end position in a read

	AlignRange(){} 
	AlignRange(int pb, int pe, int rb, int re ) {
		path_beg = pb, path_end = pe;
		read_beg = rb, read_end = re;
	}
};

typedef std::tr1::unordered_map<ReadId, AlignRange> ReadRangesMap;
typedef std::tr1::unordered_map<ReadId, bool> ReadFlagMap;

/**
 * \brief Reads placement to a path
 */
class ReadAligner
{
 private:
	GraphPath*          gpath;
    ReadRangesMap       ranges;
    ReadPlacementList   places;
	ReadFlagMap         dropped;
	GSA*                gsa;
	char**              seqs;
	PathId*             preads;
	ReadStartCount*     read_starts;
	std::string         reference;
    std::vector<size_t> depths;
	RangeCountMatrix    counts;
	int                 ltrim;
	int                 rtrim;
	int                 lstop;
	int                 rstop;
	int                 size;
	int                 nreads;
	bool                status;
	ReadAlignLog*       log;
	UsedBoundMap        used_begs;
	UsedBoundMap        used_ends;
	
	
 public:
	ReadAligner();
	ReadAligner( GraphPath *gp, GSA *g, ReadStartCount *c, char **s, int n, PathId *u, ReadAlignLog *l);
	~ReadAligner();
	void release();

	char **getReads() { return seqs; }
	ReadPlacementList *getPlacements() { return &places; }	
	size_t getSize() { return size; }
	std::string getReference() { return reference; }

	void dump( std::fstream & );
	void load( std::fstream & );

	int getLTrim() { return ltrim; }
	int getRTrim() { return rtrim; }
	std::pair<int,int> getTrims() { return std::pair<int,int>(ltrim,rtrim); }

	int getLStop() { return lstop; }
	int getRStop() { return rstop; }
	bool getStatus() { return status; }

	int getReadCount( int i, int j );

	void dropCount();	

	std::pair<int,int> getStops() { return std::pair<int,int>(lstop,rstop); }

	void getWordFrequency( WordFreqMap &word_freqs,
						   int wsize,
						   GSA *g, 
						   char **r );

	void getWordFrequencyMP( WordFreqMap &word_freqs,
							 int wsize,
							 GSA *g, 
							 char **r,
							 int cores );


	UsedBoundMap *getReadBegBoundMap() { return &used_begs; }
	UsedBoundMap *getReadEndBoundMap() { return &used_ends; }
	

 private:
	void init( GraphPath *gp, GSA *g, ReadStartCount *c, char **r, int n, PathId *u, ReadAlignLog *l);
	void align();
	void trim();
	void count();
	void update();
	void alignReads();
	void extractReads();
	void placeReads();
	void updateAlignRange( int ki, int l, BoundArray ba);
	bool placeSingleRead( std::string &rseq,
						  AlignRange  &poss,
						  ReadPlacement &place );
	bool doAlignment(AlignSummary &summary, 
					 std::string &query, 
					 std::string &sbjct );
	void updateIndel( ReadId rid,
					  int init,
					  AlignSummary &summary,
					  std::list<Mismatch> &ilist,
					  std::list<Mismatch> &dlist);
	AlignSummary compareBase(std::string &query, std::string &sbjct);
	void printPlacements();

	void computeBaseDepth();
	void adjust( int pos );

	void updateUsedCount( UsedBoundInfoLists &used_lists, bool );
	/* void updateUsedCountSP( UsedBoundInfoLists &used_lists ); */
	/* void updateUsedCountMP( UsedBoundInfoLists &used_lists ); */
};

#endif
