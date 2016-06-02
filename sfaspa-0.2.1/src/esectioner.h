/**
 * @file      esectioner.h
 * @brief     Latch region finder for two sequences
 * @details   Extension region is searched in two sequences.
 *            From a kmer match, latch region is extended.
 *            Then, kmer filter is used to find sufficient kmer match.
 *            Two different kmer sizes are used here. One for anchoring
 *            and anothr for filtering.
 * @author    Youngik Yang
 * @version   0.2
 * @date      Written in 2013
 * @date      Modified on Fri 2014-01-31 01:22:40 PM
 * @copyright J. Craig Venter Institute.
 */

#include "sectioner.h"

//typedef std::tr1::unordered_map<KmerId, std::vector<int> > KmerPossMap;

class LatchSectioner : public Sectioner
{
 private:
	int mink;               ///< minimum matching kmers
	int off_diag;           ///< allow off diagonal lengths
	int direction;          ///< latch direction
	int anchor_kmer;        ///< anchor kmer
	int filter_kmer;        ///< latch filter k-mer
	int overlap_length;     ///< minimum overlapping region
	double filter_score;    ///< latch filter score

	const KmerPosMap *qposs;      ///< query kmer positions
	const KmerPosMap *sposs;      ///< sbjct kmer positions

    const KmerArray *query_anchor_kmers; ///< anchor kmers in query
    const KmerArray *sbjct_anchor_kmers; ///< anchor kmers in sbjct

 private:
	/** Check sufficient kmer match */
	bool passKmerFilter( Section &s );

	//==============================
	// Helper functions for latching
	//==============================
	/** Locate similar region */
	bool locate();

	/** Make a latch region from given position */
	bool makeRange( Section &sect, int spos, int qpos );

 public:
	/** Constructor */
	LatchSectioner(std::string *q, 
				   std::string *s, 
				   const KmerArray *qk, 
				   const KmerArray *sk, 
				   const KmerPosMap *qp, 
				   const KmerPosMap *sp);

	//========
	// Setters
	//======== 

	void setOffDiagonal ( int l    ) { off_diag     = l; }
	void setMinSameKmers( int m    ) { mink         = m; }
	void setDirection   ( int d    ) { direction    = d; }
	void setAnchorKmer  ( int k    ) { anchor_kmer  = k; }
	void setFilterKmer  ( int k    ) { filter_kmer  = k; }
	void setFilterScore ( double s ) { filter_score = s; }

	bool find();
};
