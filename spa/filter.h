/** 
 * \file      filter.h
 * \brief     kmer filter helper.
 * \details   This filter calculates the minimum number of kmers with mismatch rate.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __KMER_FILTER_H__
#define __KMER_FILTER_H__

#include <cmath>

/** kmer filter */
namespace filter
{
	/**
	 * Minimum number of same kmers
	 * \param Integer sequence length.
	 * \param Integer size of kmer.
	 * \param Integer number of mismatch allowed.
	 * \return Integer
	 */
	inline int minSameKmerCount( int seq_len, int kmer_size, int mismatch )
	{
		return seq_len-kmer_size+1 - kmer_size*mismatch;
	}

	/**
	 * Minimum number of same kmers
	 * \param Integer sequence length.
	 * \param Integer size of kmer.
	 * \param Double  mismatch rate
	 * \return Integer
	 */
	inline int minSameKmerCount( int seq_len, int kmer_size, double mis_rate )
	{
		double mismatch = seq_len*mis_rate;
		return int( seq_len-kmer_size+1 - kmer_size*mismatch );
	}
}

#endif
