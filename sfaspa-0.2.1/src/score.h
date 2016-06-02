/** 
 * \file      score.h
 * \brief     pair-wise score
 *            
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \bug       None.
 * \warning   None.
 * \date      Modified on Tue 2013-12-17 02:36:07 PM
 * \copyright J. Craig Venter Institute.
 */

#ifndef __SCORE_H__
#define __SCORE_H__

#include <seqan/score.h>
#include <seqan/basic.h>

enum scoringMatrices { BLOSUM30, BLOSUM62, BLOSUM80 };
const int MIN = -10000000;

/**
 * \brief Score of two sequences by base composition
 */
namespace scoring 
{
	/** 
	 * Retrieve scoring matrix (BLOSUM 62:default, BLOSUM 80, BLOSUM 30)
	 */
	inline const int* getMatrix( int type ) 
	{
		using namespace seqan;
		switch(type) {
		case BLOSUM30:
			return ScoringMatrixData_<int, AminoAcid, Blosum30_>::getData();
		case BLOSUM80:
			return ScoringMatrixData_<int, AminoAcid, Blosum80_>::getData();
		default: //BLOSUM62
			return ScoringMatrixData_<int, AminoAcid, Blosum62_>::getData();
		}
	}

	/**
	 * Get score of two amino acids in scoring matrix
	 */
	inline int getScore(const seqan::AminoAcid &aa1, const seqan::AminoAcid &aa2, int mtype) 
	{
		using namespace seqan;
		const int *matrix = getMatrix( mtype );
		int index = ordValue(aa1) * ValueSize<AminoAcid>::VALUE + ordValue(aa2);
		return matrix[index];
	}

	/**
	 * Get score of two amino acids in scoring matrix
	 */
	inline int getScore(const char &ch1, const char &ch2, int mtype) 
	{
		seqan::AminoAcid aa1 = ch1, aa2 = ch2;
		return getScore(aa1, aa2, mtype);
	}

	/**
	 * Get sum of scores of two strings in same length.
	 */
	inline int sumScore(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;

		int sum = 0;
		for ( size_t i = 0; i < seq1.length(); i++ )
			sum += getScore(seq1[i], seq2[i], mtype);

		return sum;
	}

	/**
	 * Count postive scores of two strings in same length.
	 */
	inline int countPositive(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;

		int count = 0;
		for ( size_t i = 0; i < seq1.length(); i++ )
			if ( getScore(seq1[i], seq2[i], mtype) > 0 ) count++;

		return count;
	}

	/**
	 * Count postive scores of two strings in same length.
	 */
	inline int countPositive(const std::string *seq1, const std::string *seq2, size_t beg1, size_t beg2, size_t len, int mtype) 
	{
		int count = 0;
		size_t i = beg1, j = beg2;
		for ( size_t k = 0; k < len; k++ ) {
			assert( i < seq1->size() );
			assert( j < seq2->size() );
			if ( getScore((*seq1)[i], (*seq2)[j], mtype) > 0 ) count++;
			i++;   j++;
		}
		return count;
	}
	
	/**
	 * Compute postive score ratio of two strings in same length.
	 */
	inline double positiveRate(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;
		if ( seq1.length() == 0 ) return 0;

		return countPositive(seq1, seq2, mtype)/(double)seq1.length();
	}

	/**
	 * Compute postive score ratio of two strings in same length.
	 */
	inline double positiveRate(const std::string *seq1, const std::string *seq2, size_t beg1, size_t beg2, size_t len, int mtype) 
	{
		return countPositive(seq1, seq2, beg1, beg2, len, mtype)/(double)len;
	}

}

#endif
