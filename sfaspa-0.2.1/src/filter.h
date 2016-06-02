/** 
 * \file      filter.h
 * \brief     kmer filter helper.
 * \details   This filter calculates the minimum number of kmers with mismatch rate.
 * \author    Youngik Yang
 * \version   0.2
 * \date      2010-2011
 * \date      modified:  Fri 2013-12-13 04:50:26 PM
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __KMER_FILTER_H__
#define __KMER_FILTER_H__

#include <cmath>

/** 
 * \brief Kmer filter 
 */
namespace filter
{
	// L : sequence length
	// k : size of k-mer
	// M : no. of mismatch
	// X : min. no of shared k-mer
	// X >= L-k+1 - k*M
	/**
	 * Minimum number of same kmers
	 * \param seq_len    sequence length.
	 * \param kmer_size  size of kmer.
	 * \param mismatch    number of mismatch allowed.
	 * \return Integer
	 */
	inline int minSameKmerCount( int seq_len, int kmer_size, int mismatch )
	{
		return seq_len-kmer_size+1 - kmer_size*mismatch;
	}

	/**
	 * Minimum number of same kmers
	 * \param seq_len    sequence length.
	 * \param kmer_size  size of kmer.
	 * \param mis_rate   mismatch rate
	 * \return Integer
	 */
	inline int minSameKmerCount( int seq_len, int kmer_size, double mis_rate )
	{
		/* double mismatch = round(seq_len*mis_rate); */
		/* return int( seq_len-kmer_size+1 - kmer_size*mismatch ); */
		/* return round( seq_len-kmer_size+1 - kmer_size*mismatch ); */
		double mismatch = seq_len*mis_rate;
		double minkmers = seq_len-kmer_size+1 - kmer_size*mismatch;
		return minkmers;
		//return round(minkmers);
		//return floor(minkmers);
	}

	// L : sequence length
	// k : size of k-mer
	// M : no. of mismatch
	// X : min. no of shared k-mer
	// X >= L-k+1 - k*M
	// k >= (L+1-X)/(M+1)
	/** 
	 * Maximum k-mer size decision.
	 * \param seq_len sequence length.
	 * \param min_kmer minimum kmer-filter size
	 * \param mis_rate mismatch rate
	 */
	inline int getFilterKmer( int seq_len, int min_kmer, double mis_rate ) 
	{
		double mismatch  = seq_len*mis_rate;
		double min_ksize = (seq_len+1-min_kmer) / (mismatch+1);
		//return round(max_ksize);
		return min_ksize;
	}

	/* // X >= L-k+1 -k*M */
	/* // X >= L-k+1 -k*(E*L) */
	/* // X+k-1 >= L(1-k*E) */
	/* // L <= (X+k-1)/(1-K*E) */
	/* /\**  */
	/*  * Get minimum sequence length */
	/*  * \param kmer_size size of kmer */
	/*  * \param min_kmer minimun kmer match */
	/*  * \param mis_rate mismatch rate */
	/*  *\/		 */
	/* inline int getMinSequenceLength( int kmer_size, int min_kmer, double mis_rate ) */
	/* { */
	/* 	if ( kmer_size*mis_rate >= 1 ) return -1; */
			
	/* 	double length = (min_kmer+kmer_size-1) / (1-kmer_size*mis_rate); */
	/* 	return ceil(length); */
	/* } */

	/* // L-k+1 -k*M >= X */
	/* // L-k+1 -k*(E*L) >= X */
	/* // E*K*L <= L-K+1-X */
	/* // E <= (L-K+1-X)/(K*L) */
	/* inline double getErrorRate( int seq_len, int kmer_size, int min_kmer ) */
	/* { */
	/* 	return (seq_len+kmer_size+1-min_kmer)/(kmer_size*seq_len); */
	/* } */
}

#endif
