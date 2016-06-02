/**
 * \file alignsummary.h
 * Alignment summary.
 * Summary of pairwise alignment including score, length, gaps, indels, etc.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \copyright J. Craig Venter Institute.
 */

#ifndef __ALIGNSUMMARY_H__
#define __ALIGNSUMMARY_H__

#include <ostream>

/** 
 * Alignment position index.
 */
struct AlignIndex
{
	int aln_pos; ///< alignment position
	int seq_pos; ///< query position
	int ref_pos; ///< sbjct position
	AlignIndex( int a, int s, int r ) 
	{ 
		aln_pos = a; seq_pos = s; ref_pos = r; 
	}
};

typedef std::pair<AlignIndex, AlignIndex> AlignIndexPair;
typedef std::pair<int, int> intPair;
typedef std::vector<int> intVec;
typedef std::list<AlignIndex> AlignPosList;

/**
 * Alignment score summary
 */
struct AlignSummary 
{
	int score;          ///< alignment score
	int length;         ///< alignment length
	int gap;            ///< gap count
	int match;          ///< match count
	int mismatch;       ///< mismatch count
	int positive;       ///< positive count
	double posrate;     ///< positive rate
	intVec ins;         ///< insertion index
	intVec del;         ///< deletion  index
	intPair range;      ///< alignment range
	intPair outer;      ///< out boundary
	intPair lgap;       ///< leading gaps
	intPair egap;       ///< end gaps
	intPair s1se;       ///< 1st seq start & end positions
	intPair s2se;       ///< 2nd seq start & end positions
	AlignPosList ilist; ///< insertion positions
	AlignPosList dlist; ///< deletion positions

	AlignSummary() 
	{
		score = length = gap = match = mismatch = positive = 0;
		posrate = 0;
	}
	
	/**
	 * copy constructor 
	 */
	AlignSummary(const AlignSummary &source) 
	{
		score    = source.score;
		length   = source.length;
		gap      = source.gap;
		match    = source.match;
		mismatch = source.mismatch;
		positive = source.positive;
		posrate  = source.posrate;
		ins      = source.ins;
		del      = source.del; 
		range    = source.range;
		outer    = source.outer;
		lgap     = source.lgap;
		egap     = source.egap;
		s1se     = source.s1se;
		s2se     = source.s2se;
		ilist    = source.ilist;
		dlist    = source.dlist;
	}

	/**
	 * Operator overloading
	 */
	AlignSummary& operator= (const AlignSummary &source)
	{
		if (this == &source) return *this;

		score    = source.score;
		length   = source.length;
		gap      = source.gap;
		match    = source.match;
		mismatch = source.mismatch;
		positive = source.positive;
		posrate  = source.posrate;
		ins      = source.ins;
		del      = source.del; 
		range    = source.range;
		outer    = source.outer;
		lgap     = source.lgap;
		egap     = source.egap;
		s1se     = source.s1se;
		s2se     = source.s2se;
		ilist    = source.ilist;
		dlist    = source.dlist;

		return *this;
	}
	   
		void print(std::ostream &out)
	{
		out << "\nAlignment Summmary\n";
        out << "length:" << length << "\n";
        out << "match:" << match << "\n";
        out << "mismatch:" << mismatch << "\n";
        out << "positive:" << positive << "\n";
		out << "score:" << score << "\n";
        out << "%positive:" << posrate << "\n";
        out << "insertion:" << ins.size() << "\n";
        out << "deletion:" << del.size() << "\n";
        out << "sbjct range:" << s1se.first << "\t" << s1se.second << "\n";
        out << "query range:" << s2se.first << "\t" << s2se.second << "\n";
        out << "leading gap:" << lgap.first << "\t" << lgap.second << "\n";
        out << "trailing gap:" << egap.first << "\t" << egap.second << "\n";
        out << "outer range:" << outer.first << "\t" << outer.second << "\n";
        out << "inner range:" << range.first << "\t" << range.second << "\n";
		out << "\n";
	}
};

#endif
