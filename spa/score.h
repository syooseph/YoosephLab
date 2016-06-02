#ifndef __SCORE_H__
#define __SCORE_H__

#include <seqan/score.h>
#include <seqan/basic.h>

enum scoringMatrices { BLOSUM30, BLOSUM62, BLOSUM80 };
const int MIN = -10000000;

namespace scoring 
{
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

	inline int getScore(const seqan::AminoAcid &aa1, const seqan::AminoAcid &aa2, int mtype) 
	{
		using namespace seqan;
		const int *matrix = getMatrix( mtype );
		int index = ordValue(aa1) * ValueSize<AminoAcid>::VALUE + ordValue(aa2);
		return matrix[index];
	}

	inline int getScore(const char &ch1, const char &ch2, int mtype) 
	{
		seqan::AminoAcid aa1 = ch1, aa2 = ch2;
		return getScore(aa1, aa2, mtype);
	}

	inline int sumScore(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;

		int sum = 0;
		for ( size_t i = 0; i < seq1.length(); i++ )
			sum += getScore(seq1[i], seq2[i], mtype);

		return sum;
	}

	inline int countPositive(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;

		int count = 0;
		for ( size_t i = 0; i < seq1.length(); i++ )
			if ( getScore(seq1[i], seq2[i], mtype) > 0 ) count++;

		return count;
	}

	inline double positiveRate(const std::string &seq1, const std::string &seq2, int mtype) 
	{
		if ( seq1.length() != seq2.length() ) 
			return MIN;
		if ( seq1.length() == 0 ) return 0;

		return countPositive(seq1, seq2, mtype)/(double)seq1.length();
	}
}

#endif
