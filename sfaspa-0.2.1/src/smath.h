/** 
 * @file smath.h
 * @date 2013
 * @date Modified: Fri 2013-12-13 06:41:10 PM
 * @author Youngik Yang
 * @version 0.2
 * @brief Tiny math functions
 */

#ifndef __SMALL_MATH_H__
#define __SMALL_MATH_H__

#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>
#include <cassert>

/**
 * \brief Simple math functions
 */
namespace math 
{
	/** 
	 * Sum of values 
	 * \param data an array of values
	 * \param n size of the array
	 */
	double sum(const double *data, size_t n);
	
	/** 
	 * Mean of values 
	 * \param data an array of values
	 * \param n size of the array
	 */
	double mean(const double *data, size_t n); 
	
	/** 
	 * Maximum value
	 * \param data an array of values
	 * \param n size of the array
	 */
	double max(const double *data, size_t n);

	/** 
	 * Minimum value
	 * \param data an array of values
	 * \param n size of the array
	 */
	double min(const double *data, size_t n);

	/** 
	 * Median value
	 * \param data an array of values
	 * \param n size of the array
	 * \param sorted is data sorted?
	 */
	double median(double *data, size_t n, bool sorted);
	
	/** 
	 * Quantile value
	 * \param data an array of values
	 * \param n size of the array
	 * \param f quantile of interest (between 0 to 1.0)
	 * \param sorted is data sorted?
	 */
	double quantile(double *data, size_t n, double f, bool sorted);

	/** 
	 * Sort data
	 * \param data an array of values
	 * \param n size of the array
	 */
	void sort(double *data, size_t n);
}

#endif
