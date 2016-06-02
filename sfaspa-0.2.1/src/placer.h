/** 
 * @file       placer.h
 * @brief      Meta-path aligner
 * @details    Alignment of reads to meta-path (extended/clusterd).
 * @date       Thu 2014-01-02
 * @author     Youngik Yang
 * @version    0.02
 * @copyright  J. Craig Venter Institute.
 */

#ifndef __PLACER_H__
#define __PLACER_H__

#include "metapath.h"

enum invalidPath { ZERO_READS, EMPTY_SEQUENCE, TRIM_ERROR, COVERAGE_HOLE };

/**
 * \brief Read placement to meta-path.
 */
class Placer : public MetaPath
{
 private:
	size_t ct_gappy;
	size_t ct_gappy_merge;
	size_t ct_gappy_extend;

	size_t ct_realign, ct_empty_ref, ct_zero_reads;
	size_t ct_msa_empty, ct_msa_lerr, ct_msa_rerr, ct_msa_hole;
	size_t ct_con_lerr, ct_con_rerr, ct_con_hole;
	size_t ct_msa_ltrim, ct_msa_rtrim, ct_con_ltrim, ct_con_rtrim;

 private:
	void align();
	void alignAugmentPath( size_t nid, PathId pid, PathAligner &raln );
	void alignFinalMergers( size_t nid, PathId pid, PathAligner &raln );
	void alignSuperCluster( size_t nid, PathId cpivot, PathAligner &paln );

	void alignReadExtensionPath( PathId tpivot,
								 PathAligner &taln,
								 bool &gappy );
	void alignTinyExtensionPath( PathId tpivot,
								 PathAligner &taln,
								 bool &gappy );
	void alignLongExtensionPath( PathId lpivot,
								 PathAligner &laln,
								 bool &gappy ) ;

	void alignReadsJoiners( size_t nid, PathId pid, PathAligner &raln );
	void alignFinalJoiners( size_t nid, PathId pid, PathAligner &raln );
	void alignShortExtensionCluster( PathId ext_pid,
									 bool &gappy,
									 PathAligner &raln );
	void alignLongExtensionCluster( PathId ext_pid,
								 bool &gappy,
								 PathAligner &raln );

	void alignCluster( PathAligner &caln, 
					   PathId pid, 
					   int pos, 
					   AlignSummary &sum, 
					   bool & );
	void alignPath( PathAligner& raln,
					PathId pid );//, 
					/* AlignSummary &sum */
					/* ); */


	AlignSummary findAlignSummary( PathId pid, PathAlignMap &aligns);


	bool findPairRead( ReadId rid, PathId pid );

	void check();
	int findBadPaths();
	bool validPath(int i, int &why);
	void countInvalidPaths();

 public:
	void run();
};


#endif
