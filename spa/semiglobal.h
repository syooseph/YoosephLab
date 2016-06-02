/**
 * \file      semiglobal.h
 * \brief     Global alignment.
 * \details   Global pairwise alignment.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SEMI_GLOBAL_ALIGNMENT_H__
#define __SEMI_GLOBAL_ALIGNMENT_H__

#include "score.h"
#include "alignment.h"

enum ANCHORS { ANCHOR_CENTER, ANCHOR_LEFT, ANCHOR_RIGHT };

typedef seqan::Align<TString> TAlign;
typedef std::pair<int, int> IntPair;

/**
 * Global pairwise alignment 
 */
class SemiGlobalAlign : public AlignPair
{
 private:
	TAlign aln; ///< alignment type

 public:
	SemiGlobalAlign();
	SemiGlobalAlign(TString&, TString&, int);
	SemiGlobalAlign(std::string&, std::string&, int);
	SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop);
	
	TAlign getAlignment() { return aln; }
	void printGaps();
	void printAlignment(std::ostream&);
 private:
	void align(TString &, TString &, int);
	void align(std::string &, std::string &, int);
	void setSummary();
	void alignmentInnerRange();
	void alignmentOuterBoundary();
	void setSequenceRange();
	void setLeadingGaps();
	void setEndGaps();
	void setAlignmentStats();
	void setGapPositions();
	//bool adjustHead();
	void refine();
	int refineHead();
	void refineTail(int);
	void __trimIndels( int s, int e );
	void __adjustIndels( int nshift );
};


#endif

