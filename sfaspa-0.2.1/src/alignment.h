/**
 * \file alignment.h
 * \brief     Pairwise alignment.
 * \details   Base class of pairwise alignment.
 * \author    Youngik Yang
 * \version   0.2
 * \date      Written in 2010
 * \date      Modified on Thu 2013-10-24 12:54:38 PM
 * \copyright J. Craig Venter Institute.
 */
#ifndef __PAIR_ALIGNMENT_H__
#define __PAIR_ALIGNMENT_H__

#include <iostream>
#include <seqan/align.h> 
#include "alignsummary.h"

/** Alignment types */
enum ALIGN_TYPE { GLOBAL, LOCAL, GLOCAL };

/** Alignment score types */
enum SCORE_TYPE { SIMPLE, MATRIX };

/** SeqAn type helpers */
typedef seqan::String<seqan::AminoAcid> TString;
typedef seqan::Align<TString> TAlign;

/** 
 * \brief Base class for pairwise alignment
 */
class AlignPair
{
 protected:
	int score_type;        ///< Alignment score type: simple integer or scoring matrix
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
		match      = 3; 
		mismatch   = -3;
		gapopen    = -11;
		gapext     = -1; 
	}
	
	/**
	 * Destructor
	 */
	virtual ~AlignPair() {}

	/**
	 * \return Alignment summary
	 */
	AlignSummary getSummary() { return summary; }

	/**
	 * Print alignment
	 * \param out ostream object
	 */
	virtual void printAlignment(std::ostream &out) {}

	/**
	 * Set gap extenstion and open penalties
	 * \param ext gap extension penalty
	 * \param open gap open penalty
	 */
	void setGapPenalty(int ext, int open) 
	{
		gapext = ext; gapopen= open;
	}

	/**
	 * Set match/mismatch scores 
	 * \param m match score
	 * \param mm mismatch penalty
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
