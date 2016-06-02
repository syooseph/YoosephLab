/** 
 * @file       palign.h
 * @brief      Single path and its member reads placements
 * @date       Modified on Tue 2013-12-17 06:46:57 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __PATH_ALIGNER_H__
#define __PATH_ALIGNER_H__

#include "ralign.h"

/**
 * \brief A single path and placement of reads to the path
 */
class PathAligner
{
 private:
	std::string sequence;
	std::string consensus;
    ReadPlacementList places;

	int lstop,rstop;
	int ltrim,rtrim;

	bool dirty;

 private:
	bool alignReadToPath( AlignSummary &sum,
						  std::string &pstr,
						  std::string &rstr,
						  //ReadPlacement &p,
						  ReadPlacementList::iterator &it,
						  int rbeg,
						  int qbeg);
	
 public:
	PathAligner();
	PathAligner( ReadAligner &raln );
	void reset();
	std::string       getSequence()   { return sequence; }
	std::string       getConsensus()   { return consensus; }
	ReadPlacementList *getPlacements() { return &places; }
	void setSequence( std::string &s ) { sequence = s; }
	void setConsensus( std::string &c ) { consensus = c; }
	void setPlacements( ReadPlacementList &p ) { places = p; }
	void addGap( int pos );
	bool dropGaps( int pos );
	void adjust( int pos );
	void addPlacement( ReadPlacement &p ) { places.push_back(p); }
	void merge( PathAligner &other );
	void realign(char**);
	void printPlacement( std::ostream &out );
	void dump( std::fstream &out );
	void load( std::fstream &in );
	std::vector<size_t> computeBaseDepth( char **seqs );
	size_t getSize() { return places.size(); }

	void setLStop(int v) { lstop = v; }
	void setRStop(int v) { rstop = v; }
	void setLTrim(int v) { ltrim = v; }
	void setRTrim(int v) { rtrim = v; }
	void setDirty(bool v) { dirty =v; }

	int getLStop() { return lstop; }
	int getRStop() { return rstop; }
	int getLTrim() { return ltrim; }
	int getRTrim() { return rtrim; }
	bool getDirty() { return dirty; }
};


#endif
