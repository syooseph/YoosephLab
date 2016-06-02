/**
 * \file alignsummary.h
 * \brief     Pair-wise alignment summary.
 * \details   Summary of pairwise alignment including score, length, gaps, indels, etc.
 * \author    Youngik Yang
 * \version   0.2
 * \date      Written in 2010
 * \date      Modified on Fri 2013-12-13 02:15:17 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __ALIGNSUMMARY_H__
#define __ALIGNSUMMARY_H__

#include <list>
#include <vector>
#include "alignindex.h"

typedef std::pair<AlignIndex, AlignIndex> AlignIndexPair;
typedef std::pair<int, int> intPair;
typedef std::vector<int> intVec;
typedef std::list<AlignIndex> AlignPosList;

/**
 * \class AlignSummary
 * \brief Summary for pair-wise alignment
 */
class AlignSummary 
{
 public:
	int score;          ///< alignment score
	int length;         ///< alignment length
	int match;          ///< match count
	int mismatch;       ///< mismatch count
	int positive;       ///< positive count
	double posrate;     ///< positive rate
	intPair range;      ///< alignment range
	intPair outer;      ///< out boundary
	intPair lgap;       ///< leading gaps
	intPair egap;       ///< end gaps
	intPair s1se;       ///< 1st seq start & end positions
	intPair s2se;       ///< 2nd seq start & end positions
	AlignPosList ilist; ///< insertion positions
	AlignPosList dlist; ///< deletion positions

	/**
	 * Constructor
	 */
	AlignSummary();

	/**
	 * Copy constructor 
	 * \param source another AlignSummary object
	 */
	AlignSummary(const AlignSummary &source);

	/**
	 * Operator overloading
	 * \param source another AlignSummary object
	 */
	AlignSummary& operator= (const AlignSummary &source);

	/**
	 * Initalization
	 */
	void init();

	/** 
	 * Shift indel positions
	 * \param offset number of bases to shift
	 * \param ref_seq which sequence (sbjct:true or query:false) to shift 
	 */
	void shift(int offset, bool ref_seq);

	/** 
	 * Print alignment summary
	 */
	void print(std::ostream &out);

	/** 
	 * Print gaps
	 */
	void printGap(std::ostream &out);

	/**
	 * Dump a binary object to a file
	 */
	void dump( std::fstream &out ) ;

	/**
	 * Load a binary object from a file
	 */
	void load( std::fstream &in ) ;

	void self( size_t n );
};

#endif
