#ifndef __COVERAGE_H__
#define __COVERAGE_H__

#include <iostream>
#include "InvertedIndex.h"
#include "smath.h"

namespace cov
{
	inline void coverageSummary( InvertedIndex &iindex, std::ostream &out )
	{
		out << "\nCoverage Summary\n";
		InvertedIndexMap imap = iindex.getInvertedIndex();
		int n = imap.size();
		out << "\t# kmers:" << n << "\n";
		double *coverage = new double[n];

		int i = 0;
		for ( i = 0; i < n; i++ ) coverage[i] = 0;

		i = 0;
		for (InvertedIndexMap::iterator it = imap.begin(); it != imap.end(); ++it )
			coverage[i++] = it->second->size;

		math::sort(coverage, n);
		out << "\tavg:" << math::mean(coverage, n) << "\n";
		out << "\tmax:"  << math::max(coverage, n) << "\n";
		out << "\tmin:"  << math::min(coverage, n) << "\n";
		out << "\tmed:" << math::median(coverage, n, 1) << "\n";

		out << "\tPercentiles:\n";
		out << "\t";
		for ( int i = 10; i <= 90; i += 10  ) {
			out << i << ":" << math::quantile(coverage, n, i/100.0, 1) << " ";
		}
		out << "\n";
		out << "\t";
		for ( int i =91; i <= 99; i++ ) {
			out << i << ":" << math::quantile(coverage, n, i/100.0, 1) << " ";
		}
		out << "\n";
		out << "\t";
		for ( double i =99.1; i <= 100; i+=0.1 ) {
			out << i << ":" << math::quantile(coverage, n, i/100.0, 1) << " ";
		}
		out << "\n";

		delete[] coverage;
	}

}

#endif
