/** 
 * @file       progress.h
 * @brief      Show progress.
 * @date       Tue 2013-12-17 04:17:12 PM
 * @author     Youngik Yang
 * @version    0.2
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __PROGRESS_H__
#define __PROGRESS_H__

#include "timer.h"

/**
 * \brief Progress display
 */
struct Progress
{
	double ratio; ///< percent progressed
	size_t count; ///< processed count
	size_t total; ///< total size
	double stime; ///< start time

	Progress()
	{
		init();
	}

	Progress( double r, size_t c, size_t t, double s)
	{
		ratio = r;
		count = c;
		total = t;
		stime = s;
	}

	void init()
	{
		ratio = 0;
		count = 0;
		total = 0;
		stime = mytime();
	}

	void print()
	{
		fprintf( stdout, "%zu/%zu (%.0f%%) %.2f sec\n", count, total, ratio, mytime()-stime );
	}

	void showProgress()
	{
		if ( count/(double)total*100 >= ratio ) {             
			print();
			ratio += 1;
		}    
	}
	
	void reset()
	{
		init();
	}
};

#endif

	
