/** 
 * \file      eval.h
 * \brief     Evalution functions.
 * \details   This header includes entropy function.
 * \author    Youngik Yang
 * \version   0.001
 * \date      2010-2011
 * \bug       Not known.
 * \warning   None.
 * \copyright J. Craig Venter Institute.
 */

#ifndef __EVAL_H__
#define __EVAL_H__

#include <cmath>
#include <iostream>
#include <vector>
#include <map>

namespace eval
{

	
	/**
	 * Entropy calculation.
	 * \param  Count of each value.
	 * \param  Sum of all occurrence.
	 * \return Enropy
	 */
	template<typename T>
		double __calcEntropy( std::map<T, int> &vmap, int sum )
		{
			double h = 0;
			typename std::map<T, int>::iterator it;
			for ( it = vmap.begin(); it != vmap.end(); ++it ) {
				double pi = (it->second)/(double)sum;
				double lp = log2(pi);
				h -= (pi*lp);
			}
			return h;
		}
		
	/**
	 * Entropy caller for count map input.
	 */
	template<typename T>
		inline double entropy( std::map<T, int> &vmap, int sum )
		{
			return __calcEntropy(vmap, sum);
		}

	/**
	 * Utility function for computing entropy.
	 * \param  Vector  vector of values
	 * \param  Integer size of vector
	 * \return Count of each value occurrence
	 */
	template<typename T>
		std::map<T, int> __getCountMap(T *vector, int size)
	{
		std::map<T, int> vmap;
		for ( int i = 0; i < size; i++ ) {
			if ( vmap.find(vector[i]) == vmap.end() )
				vmap.insert( std::pair<T, int>(vector[i], 0) );
			vmap[vector[i]]++;
		}
		return vmap;
	}
	
	/**
	 * Entropy caller for input vector.
	 */
	template<typename T>
		inline double entropy( T* vector, int size )
		{
			std::map<T, int> vmap = __getCountMap<T>(vector, size);
			return __calcEntropy(vmap, size);
		}
}

#endif
