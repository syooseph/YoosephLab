/**
 * \file alignment.h
 * \brief     Pairwise alignment.
 * \details   Base class of pairwise alignment.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \copyright J. Craig Venter Institute.
 */
#ifndef __PAIR_ALIGNMENT_H__
#define __PAIR_ALIGNMENT_H__

#include <iostream>
#include <seqan/align.h> 
#include "alignsummary.h"

//enum ALIGN_TYPE { GLOBAL, LOCAL };
enum SCORE_TYPE { SIMPLE, MATRIX };

typedef seqan::String<seqan::AminoAcid> TString;
typedef seqan::Align<TString> TAlign;
typedef seqan::StringSet<TString, seqan::Dependent<> > TStringSet; 
typedef seqan::StringSet<TString, seqan::Dependent<> > TDepStringSet;
typedef seqan::Graph<seqan::Alignment<TDepStringSet> > TAlignGraph;

/** 
 * Base class of pairwise alignment
 */
class AlignPair
{
 protected:
	int score_type;        ///< Alignment score type: simple fixed integer or scoring matrix
	int match;             ///< Match count
	int mismatch;          ///< Mismatch count
	int gapopen;           ///< Gap open penalty
	int gapext;            ///< Gap extension penalty
	AlignSummary summary;  ///< Alignment score summary

 public:
	/** 
	 * Constructor
	 */
	AlignPair() 
	{ 
		score_type = SIMPLE;
		match    = 3; 
		mismatch = -3;
		gapopen  = -11;
		gapext   = -1; 
	}
	
	/**
	 * Destructor
	 *
	 */
	virtual ~AlignPair() {}

	/**
	 * \return Alignment summary
	 */
	AlignSummary getSummary() { return summary; }

	/**
	 * Print alignment
	 */
	virtual void printAlignment(std::ostream &) {}

	/**
	 * Set gap extenstion and open penalties
	 */
	void setGapPenalty(int ext, int open) 
	{
		gapext = ext; gapopen= open;
	}

	/**
	 * Set match/mismatch scores 
	 */
	void setMatchScore(int m, int mm ) 
	{
		match = m; mismatch = mm;
	}
	
	/**
	 * \return gap penalties.
	 */
	std::pair<int, int> getGapPenalty() { return std::pair<int, int>(gapext, gapopen); }

	/**
	 * \return match/mismatch scores
	 */
	std::pair<int, int> getMatchScore() { return std::pair<int, int>(match, mismatch); }
};

#endif
