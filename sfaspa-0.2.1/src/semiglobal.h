/**
 * \file      semiglobal.h
 * \brief     Semi-global alignment.
 * \details   Semi-global pairwise alignment.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      Updated: 2013
 * \pre       Seqan v1.4 
 * \warning   Gap handling is different between Seqan v1.3 and v1.4.
 * \copyright J. Craig Venter Institute.
 */


/* Seqan v1.4 is different from v1.3 when handle gaps 
 * 1. countGaps can be accessibly only by iterator
 * 2. beginPosition & clippedBeginPosition returns
 *    different values.
 */


#ifndef __SEMI_GLOBAL_ALIGNMENT_H__
#define __SEMI_GLOBAL_ALIGNMENT_H__

#include "score.h"
#include "alignment.h"

enum ANCHORS { ANCHOR_CENTER, ANCHOR_LEFT, ANCHOR_RIGHT };


typedef seqan::Align<TString> TAlign;

/**
 * \brief Semi-global pairwise alignment 
 */
class SemiGlobalAlign : public AlignPair
{
 private:
	TAlign aln;   ///< alignment type

	int lower;    ///< lower diagonal band
	int upper;    ///< upper diagonal band

	bool banded;  ///< banded alignment
	bool verbose; ///< verbose option

//========================
// Public member functions
//========================
 public:
	/** Default constructor */
	SemiGlobalAlign();

	/** 
	 * Constructor 
	 * \param sbjct  sbjct sequence of SeqAn type
	 * \param query  query sequence of SeqAn type
	 * \param anchor alignment region (left, right, center)
	 */
	SemiGlobalAlign(TString& sbjct, TString& query, int anchor);

	/** 
	 * Constructor 
	 * \param sbjct  sbjct sequence of STL type
	 * \param query  query sequence of STL type
	 * \param anchor alignment region (left, right, center)
	 */
	SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor);

	/** 
	 * Constructor 
	 * \param sbjct  sequence
	 * \param query  sequence
	 * \param anchor alignment region (left, right, center)
	 * \param gex    gap extension penalty
	 * \param gop    gap open penalty
	 */
	SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop);

	/** 
	 * Constructor 
	 * \param sbjct  sequence
	 * \param query  sequence
	 * \param anchor alignment region (left, right, center)
	 * \param gex    gap extension penalty
	 * \param gop    gap open penalty
	 * \param ldiag  lower diagnoal band
	 * \param udiag  upper diagnoal band
	 */
	SemiGlobalAlign(std::string &sbjct, std::string &query, int anchor, int gex, int gop, int ldiag, int udiag);

	/** Get alignment */ 
	TAlign getAlignment() { return aln; }

	/** Print all gaps */
	void printGaps(std::ostream &);

	/** Print alignment */
	void printAlignment(std::ostream&);

	/** Set verbose option */
	void setVerbose( bool v ) { verbose = v; } 

//========================
// Private member functions
//========================
 private:
	/** Align a pair of sequences */
	void align(TString &sbjct, TString &query, int anchor);

	/** Align a pair of sequences */
	void align(std::string &sbjct, std::string &query, int anchor);

	/** Set alignment summary */
	void setSummary();

	/** 
	 * Set alignment range.
	 * True alignment region by ignoring leading/trailing gaps
	 */
	void alignmentInnerRange();

	/** Set alignment regions of each sequence. */
	void setSequenceRange();

	/** Set leading gaps of each sequence */
	void setLeadingGaps();
	
	/** Set trailing gaps of each sequence */
	void setEndGaps();

	/** Set alignment statistics */
	void setAlignmentStats();

	/** Mark gap locations */
	void setGapPositions();
	
	/** Refine alignment region */
	void refine();

	/** Refine head region of alignment */
	int refineHead();

	/** Refine tail region of alignment */
	int refineTail();

	/** Trim invalid indels after refinment */
	void trimIndels(int, int);
};


#endif

