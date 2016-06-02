//  The extend alignment class can be used to continue the alignment from
//  an existing alignment.
//  
//  Sample usage:
//
//  string seed_seqA = "ACGT";  // set initial sequence
//	string seed_seqB = "GTCA"; 
//	string extend_seqA = "CGACGTAGCTAC";  // set extend sequence
//	string extend_seqB = "CGACTGACT"
//	vector<int> init_score_row(5);        // define arbitrary alignment scores 
//	init_score_row = {2, 1, 0, -4, 5};
//	vector<int> init_score_column(5);
//	init_score_column = {2, 3, -1, 4, 5};
//	EdgeScoreBundle<int> in_score;      
//	in_score.row = init_score_row;
//	in_score.row_del = init_score_row;
//	in_score.row_ins = init_score_row;
//	in_score.column = init_score_column;
//	in_score.column_del = init_score_column;
//	in_score.column_ins = init_score_column;
//	in_score.begin_row = 0;
//	in_score.begin_column = 0;
//	ScoringFunction<int> simple_score_extend; // instantiate scoring function
//	SeqAlignExtend<int> simple_align_extend(  // instantiate extend alignment
//		seed_seqA, seed_seqB, extend_seqA, extend_seqB,
//		in_score, &simple_score_extend, GLOBAL
//	);
//	simple_align_extend.Align();              // do the alignment

#include "seq_align.h"
#include "scoring_function.h"

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <tuple>

/*
	TODO: impelmenting the local and semiglobal alignment modules
*/

/*
	This class is used to extend an existing alignment.
	It requires the boundary scores of the previous alignment as input.
	
	In order to save space,
	it tries to access the DP matrix by implementing a pseudo array. 
*/

#ifndef _SEQ_ALIGN_EXTEND_H_
#define _SEQ_ALIGN_EXTEND_H_

template<typename ALIGNSCORETYPE>
class SeqAlignExtend: public SeqAlign<ALIGNSCORETYPE>	{
 public:
	SeqAlignExtend(void);
	SeqAlignExtend(
		std::string in_seed_seqA, std::string in_seed_seqB,
		std::string in_extend_seqA, std::string in_extend_seqB,
		EdgeScoreBundle<ALIGNSCORETYPE> in_edge_score,
		ScoringFunction<ALIGNSCORETYPE>* in_score_scheme,
		enum AlignMode in_mode
	);
	SeqAlignExtend(
		std::string in_seed_seqA, std::string in_seed_seqB,
		std::string in_extend_seqA, std::string in_extend_seqB,
		EdgeScoreBundle<ALIGNSCORETYPE> in_edge_score,
		ScoringFunction<ALIGNSCORETYPE>* in_score_scheme,
		const enum AlignMode in_mode,
		const unsigned int in_downward_band,
		const unsigned int in_rightward_band,
		const bool to_check_setup
	);
	~SeqAlignExtend();
	void init_basic_parameters(
	    std::string in_seed_seqA, std::string in_seed_seqB, 
	    std::string in_extend_seqA, std::string in_extend_seqB,
		  ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, enum AlignMode in_mode
  );
	
	void TraceBack(void);
	
	using SeqAlign<ALIGNSCORETYPE>::set_max_score_index;
	using SeqAlign<ALIGNSCORETYPE>::Align;
	using SeqAlign<ALIGNSCORETYPE>::AccessDPScoreMatrix;
	using SeqAlign<ALIGNSCORETYPE>::GetAlignmentScore;
	using SeqAlign<ALIGNSCORETYPE>::GetBandIntersection;
	using SeqAlign<ALIGNSCORETYPE>::GetRowScore;
	using SeqAlign<ALIGNSCORETYPE>::GetColumnScore;
	using SeqAlign<ALIGNSCORETYPE>::GetBestEdgeScore;
	using SeqAlign<ALIGNSCORETYPE>::GetBestGlobalScore;
	using SeqAlign<ALIGNSCORETYPE>::GetAlignment;
	using SeqAlign<ALIGNSCORETYPE>::SetMode;
	//void PrintSeq(void);
	
 protected:
	using SeqAlign<ALIGNSCORETYPE>::SeqA;
	using SeqAlign<ALIGNSCORETYPE>::SeqB;
	using SeqAlign<ALIGNSCORETYPE>::DPMatrixDimensionRow;
	using SeqAlign<ALIGNSCORETYPE>::DPMatrixDimensionColumn;
	using SeqAlign<ALIGNSCORETYPE>::SeqA_align;
	using SeqAlign<ALIGNSCORETYPE>::SeqB_align;
	using SeqAlign<ALIGNSCORETYPE>::Symbol_align;
	using SeqAlign<ALIGNSCORETYPE>::ScoreScheme;
	using SeqAlign<ALIGNSCORETYPE>::Mode;
	using SeqAlign<ALIGNSCORETYPE>::Banded;
	using SeqAlign<ALIGNSCORETYPE>::DownwardBand;
	using SeqAlign<ALIGNSCORETYPE>::RightwardBand;
	using SeqAlign<ALIGNSCORETYPE>::trace_back_called;
	using SeqAlign<ALIGNSCORETYPE>::final_alignment_query_begin;
	using SeqAlign<ALIGNSCORETYPE>::final_alignment_query_end;
	using SeqAlign<ALIGNSCORETYPE>::final_alignment_target_begin;
	using SeqAlign<ALIGNSCORETYPE>::final_alignment_target_end;
	using SeqAlign<ALIGNSCORETYPE>::alignment_score;
	using SeqAlign<ALIGNSCORETYPE>::Nucleotide_matching;
	
	// indicator of which parts of alignment scores are valid
	unsigned int ValidRowBegin;	// from ValidRowBegin to EdgeScoreRow.size() is valid
	unsigned int ValidColumnBegin;	// from ValidColumnBegin to EdgeScoreColumn.size() is valid
	// the edge score from the previous alignment that we try to extend
	std::string SeedSeqA;
	std::string SeedSeqB;
	std::string ExtendSeqA;
	std::string ExtendSeqB;
	std::vector<ALIGNSCORETYPE> EdgeScoreRow;
	std::vector<ALIGNSCORETYPE> EdgeScoreColumn;
	std::vector<ALIGNSCORETYPE> EdgeScoreRow_DEL;
	std::vector<ALIGNSCORETYPE> EdgeScoreColumn_DEL;
	std::vector<ALIGNSCORETYPE> EdgeScoreRow_INS;
	std::vector<ALIGNSCORETYPE> EdgeScoreColumn_INS; 
	// the scoring matrix are separated into two halfs (the upper half and lower half),
	// Define the lengths of the original sequences to be M and N, and 
	// the lengths of the sequence to be extended to be m and n.
	// The space required is thus M*n + N*m + m*n, instead of M*N + M*n + N*m + m*n.
	// M*N space is saved.
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_Up;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_Down;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_DEL_Up;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_DEL_Down;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_INS_Up;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPScoreMatrix_INS_Down;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPTraceMatrix_Up;
	std::vector< std::vector<ALIGNSCORETYPE> > pseudoDPTraceMatrix_Down;
	// dummy variable foo, out-of-range pointers point to here
	ALIGNSCORETYPE foo;
	void init_alignment(void);
	void align(void);
	inline void fill_initial_value(void);
	inline void fill_DPMatrix_row(const unsigned int row);
	using SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix_DEL;
	using SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix_INS;
	using SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix;
	using SeqAlign<ALIGNSCORETYPE>::fill_cell_DPTraceMatrix;
	using SeqAlign<ALIGNSCORETYPE>::compute_banded_region;
	using SeqAlign<ALIGNSCORETYPE>::check_band_setup;
	inline ALIGNSCORETYPE* get_cell_address(
	    const unsigned int i, const unsigned int j, const enum DPScoreMatrixName name
	);
	using SeqAlign<ALIGNSCORETYPE>::get_character;	
	using SeqAlign<ALIGNSCORETYPE>::get_global_alignment_score_best;
	using SeqAlign<ALIGNSCORETYPE>::get_band_intersection;
};

/*
template<typename ALIGNSCORETYPE> void SeqAlignExtend<ALIGNSCORETYPE>::PrintSeq(void)	{
	std::cout << SeqA << std::endl;
	std::cout << SeqB << std::endl;
	return;
}
*/

template<typename ALIGNSCORETYPE> 
SeqAlignExtend<ALIGNSCORETYPE>::SeqAlignExtend(void) : SeqAlign<ALIGNSCORETYPE>()	{
	return;
}

template<typename ALIGNSCORETYPE> 
SeqAlignExtend<ALIGNSCORETYPE>::~SeqAlignExtend(void)	{
	return;
}

template<typename ALIGNSCORETYPE> 
void SeqAlignExtend<ALIGNSCORETYPE>::init_basic_parameters(
	std::string in_seed_seqA, std::string in_seed_seqB, 
	std::string in_extend_seqA, std::string in_extend_seqB,
	ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, enum AlignMode in_mode
)	{
	assert(in_extend_seqA.length() > 0 || in_extend_seqB.length() > 0);
	
	SeedSeqA = in_seed_seqA;
	SeedSeqB = in_seed_seqB;
	ExtendSeqA = in_extend_seqA;
	ExtendSeqB = in_extend_seqB;
	
	SeqA = SeedSeqA + ExtendSeqA;
	SeqB = SeedSeqB + ExtendSeqB;
	ScoreScheme = in_score_scheme;
	Mode = in_mode;
	return;
}


template<typename ALIGNSCORETYPE> SeqAlignExtend<ALIGNSCORETYPE>::SeqAlignExtend(
		std::string in_seed_seqA, std::string in_seed_seqB,
		std::string in_extend_seqA, std::string in_extend_seqB,
		EdgeScoreBundle<ALIGNSCORETYPE> in_edge_score,
		ScoringFunction<ALIGNSCORETYPE>* in_score_scheme,
		enum AlignMode in_mode
)	: SeqAlign<ALIGNSCORETYPE>() {
	
	assert(in_extend_seqA.length() > 0 || in_extend_seqB.length() > 0);
	assert(in_seed_seqA.length() + 1 == in_edge_score.column.size());
	assert(in_seed_seqB.length() + 1 == in_edge_score.row.size());
	assert(in_edge_score.row.size() == in_edge_score.row_del.size());
	assert(in_edge_score.row.size() == in_edge_score.row_ins.size());
	assert(in_edge_score.column.size() == in_edge_score.column_del.size());
	assert(in_edge_score.column.size() == in_edge_score.column_ins.size());
	assert(*(in_edge_score.row.end() - 1) == *(in_edge_score.column.end() - 1));
	assert(*(in_edge_score.row_del.end() - 1) == *(in_edge_score.column_del.end() - 1));
	assert(*(in_edge_score.row_ins.end() - 1) == *(in_edge_score.column_ins.end() - 1));
	assert(in_edge_score.begin_row >= 0 && in_edge_score.begin_row <= in_seed_seqB.length());
	assert(in_edge_score.begin_column >= 0 && in_edge_score.begin_column <= in_seed_seqA.length());
	if(in_edge_score.begin_row >= in_seed_seqB.length())	{
		assert(in_extend_seqB.length() > 0);
	}
	if(in_edge_score.begin_column >= in_seed_seqA.length())	{
		assert(in_extend_seqA.length() > 0);
	}
	
	init_basic_parameters(in_seed_seqA, in_seed_seqB, in_extend_seqA, in_extend_seqB, in_score_scheme, in_mode);
	DPMatrixDimensionRow = in_seed_seqA.length() + in_extend_seqA.length() + 1;
	DPMatrixDimensionColumn = in_seed_seqB.length() + in_extend_seqB.length() + 1;
	
	ValidRowBegin = in_edge_score.begin_row;
	ValidColumnBegin = in_edge_score.begin_column;
	
	EdgeScoreRow = in_edge_score.row;
	EdgeScoreColumn = in_edge_score.column;
	EdgeScoreRow_DEL = in_edge_score.row_del;
	EdgeScoreColumn_DEL = in_edge_score.column_del;
	EdgeScoreRow_INS = in_edge_score.row_ins;
	EdgeScoreColumn_INS = in_edge_score.column_ins;
	
	Banded = false;
	
	return;
}


template<typename ALIGNSCORETYPE> SeqAlignExtend<ALIGNSCORETYPE>::SeqAlignExtend(
		std::string in_seed_seqA, std::string in_seed_seqB,
		std::string in_extend_seqA, std::string in_extend_seqB,
	  EdgeScoreBundle<ALIGNSCORETYPE> in_edge_score,
		ScoringFunction<ALIGNSCORETYPE>* in_score_scheme,
		const enum AlignMode in_mode,
		const unsigned int in_downward_band, const unsigned int in_rightward_band,
		const bool to_check_setup
)	: SeqAlign<ALIGNSCORETYPE>() {
	
	assert((in_extend_seqA.length() > 0)|| (in_extend_seqB.length() > 0));
	assert(in_seed_seqA.length() + 1 == in_edge_score.column.size());
	assert(in_seed_seqB.length() + 1 == in_edge_score.row.size());
	assert(in_edge_score.row.size() == in_edge_score.row_del.size());
	assert(in_edge_score.row.size() == in_edge_score.row_ins.size());
	assert(in_edge_score.column.size() == in_edge_score.column_del.size());
	assert(in_edge_score.column.size() == in_edge_score.column_ins.size());
	assert(*(in_edge_score.row.end() - 1) == *(in_edge_score.column.end() - 1));
	assert(*(in_edge_score.row_del.end() - 1) == *(in_edge_score.column_del.end() - 1));
	assert(*(in_edge_score.row_ins.end() - 1) == *(in_edge_score.column_ins.end() - 1));
	assert(in_edge_score.begin_row >= 0 && in_edge_score.begin_row <= in_seed_seqB.length());
	assert(in_edge_score.begin_column >= 0 && in_edge_score.begin_column <= in_seed_seqA.length());
	if(in_edge_score.begin_row >= in_seed_seqB.length())	{
		assert(in_extend_seqB.length() > 0);
	}
	if(in_edge_score.begin_column >= in_seed_seqA.length())	{
		assert(in_extend_seqA.length() > 0);
	}
	assert(in_downward_band > 0);
	assert(in_rightward_band > 0);
	
	init_basic_parameters(in_seed_seqA, in_seed_seqB, in_extend_seqA, in_extend_seqB, in_score_scheme, in_mode);
	DPMatrixDimensionRow = in_seed_seqA.length() + in_extend_seqA.length() + 1;
	DPMatrixDimensionColumn = in_seed_seqB.length() + in_extend_seqB.length() + 1;
	
	ValidRowBegin = in_edge_score.begin_row;
	ValidColumnBegin = in_edge_score.begin_column;
	
	EdgeScoreRow = in_edge_score.row;
	EdgeScoreColumn = in_edge_score.column;
	EdgeScoreRow_DEL = in_edge_score.row_del;
	EdgeScoreColumn_DEL = in_edge_score.column_del;
	EdgeScoreRow_INS = in_edge_score.row_ins;
	EdgeScoreColumn_INS = in_edge_score.column_ins;
	
	Banded = true;
	DownwardBand = in_downward_band;
	RightwardBand = in_rightward_band;
	// check banded region setup
	// checks the setup of the bands
	if(to_check_setup)	{
		int k = _min<int>(in_downward_band, in_rightward_band);
		check_band_setup(k);
	}
	return;
}

template<typename ALIGNSCORETYPE> void SeqAlignExtend<ALIGNSCORETYPE>::init_alignment(void)	{
	// the purpose is to construct a psedo matrix with size (k + m) * (l + n)
	// where k = EdgeScoreColumn.size() - ValidColumnBegin + 1, 
	// and l = EdgeScoreRow.size() - ValidRowBegin + 1
	
	
	unsigned int i;
	
	int up_dim_row = EdgeScoreColumn.size() - ValidColumnBegin;
	
	pseudoDPScoreMatrix_Up.resize(up_dim_row);
	pseudoDPScoreMatrix_DEL_Up.resize(up_dim_row);
	pseudoDPScoreMatrix_INS_Up.resize(up_dim_row);
	pseudoDPTraceMatrix_Up.resize(up_dim_row);
	for(i = 0; i < pseudoDPScoreMatrix_Up.size(); ++ i)	{
		pseudoDPScoreMatrix_Up[i].resize(ExtendSeqB.length() + 1, -_MAX_VALUE);
		pseudoDPScoreMatrix_DEL_Up[i].resize(ExtendSeqB.length() + 1, -_MAX_VALUE);
		pseudoDPScoreMatrix_INS_Up[i].resize(ExtendSeqB.length() + 1, -_MAX_VALUE);
		pseudoDPTraceMatrix_Up[i].resize(ExtendSeqB.length() + 1);
	}
	
	int down_dim_column = SeqB.length() + 1;
	
	pseudoDPScoreMatrix_Down.resize(ExtendSeqA.length() + 1);
	pseudoDPScoreMatrix_DEL_Down.resize(ExtendSeqA.length() + 1);
	pseudoDPScoreMatrix_INS_Down.resize(ExtendSeqA.length() + 1);
	pseudoDPTraceMatrix_Down.resize(ExtendSeqA.length() + 1);
	for(i = 0; i < pseudoDPScoreMatrix_Down.size(); ++ i)	{
		pseudoDPScoreMatrix_Down[i].resize(down_dim_column, -_MAX_VALUE);
		pseudoDPScoreMatrix_DEL_Down[i].resize(down_dim_column, -_MAX_VALUE);
		pseudoDPScoreMatrix_INS_Down[i].resize(down_dim_column, -_MAX_VALUE);
		pseudoDPTraceMatrix_Down[i].resize(down_dim_column);
	}
	
	return;
}


template<typename ALIGNSCORETYPE> 
inline ALIGNSCORETYPE* SeqAlignExtend<ALIGNSCORETYPE>::get_cell_address(
    const unsigned int i, const unsigned int j, const enum DPScoreMatrixName name
)	{

	assert(i >= 0 && j >= 0); 
	assert(i < DPMatrixDimensionRow);
	assert(j < DPMatrixDimensionColumn);
	// computes the effective address of the pseudo matrix
	
	//std::cout << "seq_align_extend::get_cell_address:  " << i << " " << j << std::endl;
	bool is_Up = true;
	int effective_i = -1, effective_j = -1;
	if(i < ValidColumnBegin || j < ValidRowBegin || 
		(i < SeedSeqA.length() && j < SeedSeqB.length()))	{
		foo = -_MAX_VALUE;
		return &foo;
	}	else if(i < SeedSeqA.length())	{
		is_Up = true;
		effective_i = i - ValidColumnBegin;
		effective_j = j - SeedSeqB.length();
	}	else	{
		is_Up = false;
		effective_i = i - SeedSeqA.length();
		effective_j = j - ValidRowBegin;
	}
  assert(effective_i >= 0 && effective_j >= 0);
	
	if(name == FULL)	{
		return (is_Up ? &pseudoDPScoreMatrix_Up[effective_i][effective_j] : &pseudoDPScoreMatrix_Down[effective_i][effective_j]);
	}	else if(name == DEL)	{
		return (is_Up ? &pseudoDPScoreMatrix_DEL_Up[effective_i][effective_j] : &pseudoDPScoreMatrix_DEL_Down[effective_i][effective_j]);
	}	else if(name == INS)	{
		return (is_Up ? &pseudoDPScoreMatrix_INS_Up[effective_i][effective_j] : &pseudoDPScoreMatrix_INS_Down[effective_i][effective_j]);
	}	else if(name == TRACE)	{
		return (ALIGNSCORETYPE*) 
			(is_Up ? &pseudoDPTraceMatrix_Up[effective_i][effective_j] : &pseudoDPTraceMatrix_Down[effective_i][effective_j]);
	}
	return NULL;
}

template<typename ALIGNSCORETYPE> inline void SeqAlignExtend<ALIGNSCORETYPE>::fill_initial_value(void)	{
	
	//std::cout << "Begin of SeqAlignExtend::fill_initial_value" << std::endl;
	
	unsigned int i;
	// copy the originial edge scores
	for(i = ValidRowBegin; i <= SeedSeqB.length(); ++ i)	{
		*get_cell_address(SeedSeqA.length(), i, FULL) = EdgeScoreRow[i];
		*get_cell_address(SeedSeqA.length(), i, DEL) = EdgeScoreRow_DEL[i];
		*get_cell_address(SeedSeqA.length(), i, INS) = EdgeScoreRow_INS[i];
		//std::cout << "initialize score: " << AccessDPScoreMatrix(SeedSeqA.length(), i, FULL) << std::endl;
	}
	for(i = ValidColumnBegin; i <= SeedSeqA.length(); ++ i)	{
		*get_cell_address(i, SeedSeqB.length(), FULL) = EdgeScoreColumn[i];
		*get_cell_address(i, SeedSeqB.length(), DEL) = EdgeScoreColumn_DEL[i];
		*get_cell_address(i, SeedSeqB.length(), INS) = EdgeScoreColumn_INS[i];
		//std::cout << "initialize score: " << AccessDPScoreMatrix(i, SeedSeqB.length(), FULL) << std::endl;
	}
	
	
	// fills the boundary scores
	// determine the boundary to be filled
	unsigned int down_bound = 
	    (Banded && DownwardBand - 1 < DPMatrixDimensionRow) ? DownwardBand - 1 : DPMatrixDimensionRow - 1;
	unsigned int right_bound = 
	    (Banded && RightwardBand - 1 < DPMatrixDimensionColumn) ? RightwardBand - 1 : DPMatrixDimensionColumn - 1;
	for(i = SeedSeqA.length() + 1; i <= down_bound; ++ i)	{
		// fill the left-most column, fill FULL and DEL only
		if(i == SeedSeqA.length() + 1)	{
			*get_cell_address(i, 0, FULL) = _max<ALIGNSCORETYPE>(
				*get_cell_address(i - 1, 0, FULL) + ScoreScheme->GetGapOpen() + ScoreScheme->GetGapExtend(),
				*get_cell_address(i - 1, 0, INS) + ScoreScheme->GetGapExtend()
			);
		}	else	{
			*get_cell_address(i, 0, FULL) = *get_cell_address(i - 1, 0, FULL) + ScoreScheme->GetGapExtend();	
		}
		*get_cell_address(i, 0, DEL) = *get_cell_address(i, 0, FULL) + ScoreScheme->GetGapOpen();
		*get_cell_address(i, 0, INS) = *get_cell_address(i, 0, FULL);
		*(enum AlignPath*)get_cell_address(i, 0, TRACE) = UP;
	}
	
	for(i = SeedSeqB.length() + 1; i <= right_bound; ++ i)	{
		// fill the up-most row, fill FULL and INS only
		if(i == SeedSeqB.length() + 1)	{
			*get_cell_address(0, i, FULL) = _max<ALIGNSCORETYPE>(
				*get_cell_address(0, i - 1, FULL) + ScoreScheme->GetGapOpen() + ScoreScheme->GetGapExtend(),
				*get_cell_address(0, i - 1, DEL) + ScoreScheme->GetGapExtend()
			);
		}	else	{
			*get_cell_address(0, i, FULL) = *get_cell_address(0, i - 1, FULL) + ScoreScheme->GetGapExtend();	
		}
		*get_cell_address(0, i, INS) = *get_cell_address(0, i, FULL) + ScoreScheme->GetGapOpen();
		*get_cell_address(0, i, DEL) = *get_cell_address(0, i, FULL);
		*(enum AlignPath*)get_cell_address(0, i, TRACE) = LEFT;
	}
	
	//std::cout << "SeqAlignExtend::fill_initial_value finish computing boundary values" << std::endl;
	
	// filling band-specific values
	if(Banded)	{
		//std::cout << "SeqAlignExtend::fill_initial_value" << std::endl;
		//std::cout << DownwardBand << "	" << RightwardBand << std::endl;
		//std::cout << DPMatrixDimensionRow << "	" << DPMatrixDimensionColumn << std::endl;
		if(DownwardBand < DPMatrixDimensionRow)	{
			*get_cell_address(DownwardBand, 0, FULL) = *get_cell_address(DownwardBand, 0, DEL) = -_MAX_VALUE;
		}
		if(RightwardBand < DPMatrixDimensionColumn)	{
			*get_cell_address(0, RightwardBand, FULL) = *get_cell_address(0, RightwardBand, INS) = -_MAX_VALUE;
		}
		//std::cout << "SeqAlignExtend::fill_initial_value boundary negative values filled" << std::endl;
		for(i = 1; i < DPMatrixDimensionRow; ++ i)	{
			//std::cout << "== " << i << " ==" << std::endl;
			std::pair <unsigned int, unsigned int> bound = compute_banded_region(i);
			//std::cout << bound.first << "	" << bound.second << std::endl;
			if(bound.first - 1 > 0 && bound.first - 1 < DPMatrixDimensionColumn)	{
			//	std::cout << "bound to fill****" << i << "\t" << bound.first - 1 << std::endl; 
				*get_cell_address(i, bound.first - 1, FULL) = *get_cell_address(i, bound.first - 1, DEL) = -_MAX_VALUE;
			//	std::cout << "finished" << std::endl;
			}
			//std::cout << "a" << std::endl;
			if(bound.second + 1 > 0 && bound.second + 1 < DPMatrixDimensionColumn)	{
				*get_cell_address(i, bound.second + 1, FULL) = *get_cell_address(i, bound.second + 1, INS) = -_MAX_VALUE;
			}
			//std::cout << "b" << std::endl;
		}
	}
	
	//std::cout << "End of SeqAlignExtend::fill_initial_value" << std::endl;
	return;
}

template<typename ALIGNSCORETYPE> 
inline void SeqAlignExtend<ALIGNSCORETYPE>::fill_DPMatrix_row(const unsigned int row)	{
	assert(row > 0 && row < DPMatrixDimensionRow);
	// check if banded alignment is needed
	unsigned int left_bound, right_bound = DPMatrixDimensionColumn - 1;
	if(row < SeedSeqA.length() + 1)	{
		left_bound = SeedSeqB.length() + 1;
	}	else	{
		left_bound = ValidRowBegin + 1;
	}
	if(Banded)	{
		std::pair<unsigned int, unsigned int> bound = compute_banded_region(row);
		left_bound = bound.first > left_bound ? bound.first : left_bound;
		right_bound = bound.second < right_bound ? bound.second : right_bound;
	}
	unsigned int i;
	//std::cout << "left and rignt bound: " << left_bound << "\t" << right_bound << std::endl;
	for(i = left_bound; i <= right_bound; ++ i)	{
		fill_cell_DPScoreMatrix_DEL(row, i);
		fill_cell_DPScoreMatrix_INS(row, i);
		fill_cell_DPScoreMatrix(row, i);
		fill_cell_DPTraceMatrix(row, i);
		//std::cout << "Alignment Score: " << AccessDPScoreMatrix(row, i, FULL) << std::endl;
	}
	return;
}	


template<typename ALIGNSCORETYPE> void SeqAlignExtend<ALIGNSCORETYPE>::align()	{
	assert(SeqA.length() > 0 && SeqB.length() > 0);
	//std::cout << "before calling fill_initial_value" << std::endl;
	fill_initial_value();
	//std::cout << "SeqAlignExtend::fill_initial_value finished" << std::endl;
	unsigned int i;
	for(i = 1; i < DPMatrixDimensionRow; ++ i)	{
	  //std::cout << i << std::endl;
		fill_DPMatrix_row(i);
	}
	//std::cout << "Alignment finished" << std::endl;
	return;
}

template<typename ALIGNSCORETYPE> 
void SeqAlignExtend<ALIGNSCORETYPE>::TraceBack(void)	{
	trace_back_called = true;
	if(Banded)	{
		std::vector<unsigned int> intersection(4);
		intersection = get_band_intersection();
		//std::cout << "!!  " << SeedSeqA << " " << SeedSeqB << std::endl;
		//std::cout << "!!  " << ExtendSeqA << " " << ExtendSeqB << std::endl;
		//std::cout << "!!  " << intersection[0] << "  " << intersection[2] << std::endl;
		//std::cout << "!!  " << DPMatrixDimensionColumn << "  " << DPMatrixDimensionRow << std::endl;
		//assert(intersection[0] < DPMatrixDimensionColumn);
		//assert(intersection[2] < DPMatrixDimensionRow);
	}
	int idx_i, idx_j;
	set_max_score_index(idx_i, idx_j);
	//std::cout << "Check dimension:  " << DPMatrixDimensionRow << "  " << DPMatrixDimensionColumn << " " << idx_i << " " << idx_j << std::endl;  
	final_alignment_query_end = final_alignment_query_begin = idx_i - 1;
	final_alignment_target_end = final_alignment_target_begin = idx_j - 1;
	alignment_score = *get_cell_address(idx_i, idx_j, FULL);
	//std::cout << "TraceBack:: indices: " << idx_i << " " << idx_j << std::endl;
	while((idx_i > (int) SeedSeqA.length() || idx_j > (int) SeedSeqB.length()) &&
	    *(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) != TERMINAL)	{
		if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == DIAGONAL)	{
			-- idx_i;
			-- idx_j;
			final_alignment_query_begin = idx_i;
	    final_alignment_target_begin = idx_j;
			SeqA_align.resize(SeqA_align.size() + 1, SeqA[idx_i]);
			SeqB_align.resize(SeqB_align.size() + 1, SeqB[idx_j]);
			Nucleotide_matching[idx_i] = idx_j;
			char sym = (*get_character(SeqA, idx_i) == *get_character(SeqB, idx_j) ? *get_character(SeqA, idx_i) : ' ');
			if(ScoreScheme->CheckMatchScore(*get_character(SeqA, idx_i), *get_character(SeqB, idx_j)) > 0 &&
			   sym == ' ')  {
			  sym = '+'; 
			}
			Symbol_align.resize(Symbol_align.size() + 1, sym);
		}	else if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == LEFT)	{
			-- idx_j;
			final_alignment_target_begin = idx_j;
			SeqA_align.resize(SeqA_align.size() + 1, '-');
			SeqB_align.resize(SeqB_align.size() + 1, SeqB[idx_j]);
			Symbol_align.resize(Symbol_align.size() + 1, ' ');
		}	else if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == UP)	{
			-- idx_i;
			final_alignment_query_begin = idx_i;
			SeqA_align.resize(SeqA_align.size() + 1, SeqA[idx_i]);
			SeqB_align.resize(SeqB_align.size() + 1, '-');
			Nucleotide_matching[idx_i] = -1;
			Symbol_align.resize(Symbol_align.size() + 1, ' ');
		}
	}
	SeqA_align = std::string(SeqA_align.rbegin(), SeqA_align.rend());
	SeqB_align = std::string(SeqB_align.rbegin(), SeqB_align.rend());
	Symbol_align = std::string(Symbol_align.rbegin(), Symbol_align.rend());
	
	return;
}


#endif
