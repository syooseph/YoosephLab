/**
 * \file      galignment.h
 * \brief     Global alignment.
 * \details   Global pairwise alignment.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \copyright J. Craig Venter Institute.
 */

#ifndef __GLOBAL_PAIR_ALIGNMENT_H__
#define __GLOBAL_PAIR_ALIGNMENT_H__

#include "score.h"
#include "alignment.h"

typedef seqan::Align<TString> TAlign;

/**
 * Global pairwise alignment 
 */
class GlobalAlignPair : public AlignPair
{
 private:
	TAlign aln; ///< alignment type

 public:
	GlobalAlignPair();
	GlobalAlignPair(TString&, TString&);
	GlobalAlignPair(std::string&, std::string&);

	TAlign getAlignment() { return aln; }
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

