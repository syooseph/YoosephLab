#include "sectioner.h"

//typedef std::multimap<int,int> PosPairMap;

class MergeSectioner : public Sectioner
{
 private:
	int        extra;       ///< extra left/right bases to make sub-sbjct wider
	int        kmer_size;   ///< k-mer size
	const KmerArray  *query_kmers; ///< query k-mers
	const KmerArray  *sbjct_kmers; ///< sbjct k-mers
	const KmerPosMap *qposs;
	const KmerPosMap *sposs;
		
	PosPairMap poss_pair;   ///< kmer match positions between two sequences

	/* bool   aln_flag; */
	/* bool   align_success; */
	/* double merge_score; */
	/* double align_score; */

 private:
	//================================
	// Helper functions for clustering
	//================================
	/** Make kmer position pair map */
	void makePairMap();

	/** Locate similar region */
	void locate();

	/** Extend region */
	void expand();


 public:
	/** Constructor */
	MergeSectioner(std::string *q, 
				   std::string *s, 
				   const KmerArray *qk, 
				   const KmerArray *sk,
				   const KmerPosMap *qp, 
				   const KmerPosMap *sp);

	/** Set extra base */
	void setExtraBases  ( int s ) { extra = s; }
	
	/** Set filter kmer size */
	void setKmerSize( int k ) { kmer_size = k; }

	/** Find region */
	bool find();

	void pad();
	/* bool aligned() { return align_success; } */

	/* double getAlignScore() { return align_score; } */
};
