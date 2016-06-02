#ifndef __SMALL_MATH_H__
#define __SMALL_MATH_H__

//#include <gsl/gsl_statistics.h>
//#include <gsl/gsl_sort.h>
#include <numeric>
#include <functional>
#include <algorithm>
#include <vector>
//#include "debruijn.h"

namespace math 
{
	inline double sum(const double *data, size_t n) 
	{
		return std::accumulate(data, data+n, 0);
	}
	/* inline double sum(const double *data, size_t n)  */
	/* { */
	/* 	double s = 0; */
	/* 	for ( size_t i = 0; i < n; i++ ) s += data[i]; */
	/* 	return s; */
	/* } */
	
	inline double mean(const double *data, size_t n) 
	{
		return sum(data,n)/n;
	}


	/* inline double mean(const double *data, size_t n)  */
	/* { */
	/* 	return gsl_stats_mean(data, 1, n); */
	/* } */

	/* inline double variance(const double *data, size_t n)  */
	/* { */
	/* 	return gsl_stats_variance(data, 1, n); */
	/* } */

	/* inline double sd(const double *data, size_t n)  */
	/* { */
	/* 	return gsl_stats_sd (data, 1, n); */
	/* } */
	
	inline double max(const double *data, size_t n)  
	{
		return *std::max_element(data, data+n);
	}
	/* inline double max(const double *data, size_t n)  */
	/* { */
	/* 	return gsl_stats_max(data, 1, n); */
	/* } */
	
	inline double min(const double *data, size_t n)  
	{
		return *std::min_element(data, data+n);
	}
	/* inline double min(const double *data, size_t n)  */
	/* {		 */
	/* 	return gsl_stats_min(data, 1, n); */
	/* } */
	
	inline double median(double *data, size_t n, bool sorted)  
	{
		if (!sorted) std::sort(data, data+n);
		return *(data+n/2);
	}
	/* inline double median(double *data, size_t n, bool sorted)  */
	/* { */
	/* 	if (!sorted) gsl_sort(data, 1, n); */
	/* 	return gsl_stats_median_from_sorted_data (data, 1, n); */
	/* } */
	
	inline double quantile(double *data, size_t n, double f, bool sorted)  
	{
		//using namespace std;
		if (!sorted) std::sort(data, data+n);
		//return std::nth_element( data, data+int(n*f), data+n );
		return *(data+int(n*f));
	}

	/* inline double quantile(double *data, size_t n, double f, bool sorted)  */
	/* { */
	/* 	if (!sorted) gsl_sort(data, 1, n); */
	/* 	return gsl_stats_quantile_from_sorted_data (data, 1, n, f); */
	/* } */
	
	inline void sort(double *data, size_t n)  
	{
		std::sort(data, data+n);
	}
	/* inline void sort(double *data, size_t n)  */
	/* { */
	/* 	gsl_sort(data, 1, n); */
	/* } */
}

#endif
