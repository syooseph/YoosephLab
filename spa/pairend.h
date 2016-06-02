#ifndef __PAIR_END_H__
#define __PAIR_END_H__

//#include "math.h"

namespace pe 
{
	/*
	  Illumina paired-end
	  ooooooooooooo     oooooooooooooo
	     ----                 ----
		     |----|XXXXX|----|
			   d1    m    d2
			 
			 _____insert_____
	 */
	int distance( std::string &lstr,
				  std::string &rstr,
				  BitString *bstrs,
				  ReadId *lreads,
				  ReadId *rreads,
				  int *linits,
				  int *rinits,
				  size_t npair,
				  int insert,
				  int platform )
	{
		if ( platform != ILLUMINA )  {
			std::cerr << "Not supporting platform\n";
			exit(1);
		}

/* 		double min_ins = insert - 2*stddev; */
/* 		double max_ins = insert + 2*stddev; */

		double *dists = new double[npair];
		for ( size_t i = 0; i < npair; i++ ) {
			std::string rstr = bstrs[lreads[i]].toString();
			int d1 = lstr.size() - ( linits[i] + rstr.length() );
			int d2 = rinits[i];
			
/* 			double min_dist = min_ins - d1 - d2; */
/* 			double max_dist = max_int - d1 - d2; */
			double avg_dist = insert  - d1 - d2;
			//std::cerr << i << "\t" << min_dist << "\t" << max_dist << "\t" << avg_dist << "\n";

			dists[i] = avg_dist;
		}
		
		double avg = math::mean(dists, npair);
		return (int)avg;
	}
}

#endif
