//==============================================================================
//Tue 2011-06-14 12:35:24 PM
//==============================================================================
#ifndef __LOCAL_PAIR_ALIGNMENT_H__
#define __LOCAL_PAIR_ALIGNMENT_H__

//#include <seqan/align.h> 
#include "score.h"
#include "alignment.h"


class LocalAlignPair  : public AlignPair
{
 private:
	TAlignGraph alignG;

 public:
	LocalAlignPair();
	LocalAlignPair(TString&, TString&);
	LocalAlignPair(std::string&, std::string&);
	LocalAlignPair(std::string &sbjct, std::string &query, int gex, int gop);
	
	TAlignGraph getAlignment() { return alignG; }
	void printGaps();
	void printAlignment(std::ostream&);
 private:
	void align(TString &, TString &);
	void align(std::string &, std::string &);
	void setSummary();
	void alignmentInnerRange();
	void alignmentOuterBoundary();
	void setSequenceRange();
	void setLeadingGaps();
	void setEndGaps();
	void setAlignmentStats();
	void setGapPositions();
};


#endif

