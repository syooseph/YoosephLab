#include "scoring_function.h"

#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <list>
#include <tuple>
#include <unordered_map>

/*
	TODO: impelmenting the local and semiglobal alignment modules
*/

#ifndef _SEQ_ALIGN_H_
#define _SEQ_ALIGN_H_

#ifndef _MAX_VALUE
#define _MAX_VALUE 30000
#endif


template <typename NUMTYPE>
inline NUMTYPE _max(NUMTYPE num_1, NUMTYPE num_2)	{
	return (num_1 > num_2 ? num_1 : num_2);
}

template<typename NUMTYPE>
inline NUMTYPE _min(NUMTYPE num_1, NUMTYPE num_2)	{
	return (num_1 > num_2 ? num_2 : num_1);
}

template<typename SCORETYPE> 
struct EdgeScoreBundle  {
  std::vector<SCORETYPE> row;
  std::vector<SCORETYPE> column;
  std::vector<SCORETYPE> row_del;
  std::vector<SCORETYPE> column_del;
  std::vector<SCORETYPE> row_ins;
  std::vector<SCORETYPE> column_ins;
  unsigned int begin_row;
  unsigned int begin_column;
};

enum AlignMode {GLOBAL, LOCAL, SEMIGLOBAL, SEMILOCAL, HAMMING};
enum AlignPath {DIAGONAL, LEFT, UP, TERMINAL};
enum DPScoreMatrixName {FULL, DEL, INS, TRACE};

template <typename ALIGNSCORETYPE>
class SeqAlign	{
 public:	
	SeqAlign();
	SeqAlign(
	    std::string &in_seqA, std::string &in_seqB, 
	    ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, 
	    enum AlignMode in_mode
	);
	SeqAlign(
	    std::string &in_seqA, std::string &in_seqB, 
	    ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, 
	    enum AlignMode in_mode, 
		  unsigned int in_downward_band, unsigned int in_rightward_band, 
		  bool to_check_setup
  );
	~SeqAlign();
	inline void Align();	
	inline void PrintAlignment();
	inline ALIGNSCORETYPE GetAlignmentScore();
	// Interface for accessing the score in the DPScoreMatrix
	inline ALIGNSCORETYPE AccessDPScoreMatrix(
	    const unsigned int i, const unsigned int j, const enum DPScoreMatrixName name
	)  {
  	assert(i >= 0 && j >= 0);
  	assert(i < DPMatrixDimensionRow);
  	assert(j < DPMatrixDimensionColumn);
	  //std::cout << "Within AccessDPScoreMatrix:	" << i << "	" << j << std::endl;
	  if(name == TRACE)	{
	  	return (ALIGNSCORETYPE) *get_cell_address(i, j, name);
	  }	else	{
	  	return *get_cell_address(i, j, name);
	  }
  }
	inline std::vector<unsigned int> GetBandIntersection(void);
	virtual inline std::vector<ALIGNSCORETYPE> GetRowScore(unsigned int row, enum DPScoreMatrixName name);
	virtual inline std::vector<ALIGNSCORETYPE> GetColumnScore(unsigned int column, enum DPScoreMatrixName name);
	virtual inline ALIGNSCORETYPE GetBestEdgeScore();
	virtual inline ALIGNSCORETYPE GetBestGlobalScore();
	virtual inline ALIGNSCORETYPE GetBestScore();
	virtual inline void SetMode(AlignMode in_mode);
	ALIGNSCORETYPE ComputeHammingDistance(void);
	void FillEdgeScore(EdgeScoreBundle<ALIGNSCORETYPE>& edge_score);
	virtual void TraceBack(void);
	virtual std::vector<int> GetAlignmentRegion(void);
	virtual void GetAlignment(std::unordered_map<int, int>& nuc_match, std::string& seqA_align, std::string& seqB_align, std::string& symbol_align);
	void ListAllHSPPositions(
      ALIGNSCORETYPE score_cutoff, 
      std::list<std::vector<int> >& hsp_range
  );
 protected:
	// define the sequences to be aligned
	std::string SeqA;
	std::string SeqB;
	ALIGNSCORETYPE alignment_score;
	ALIGNSCORETYPE best_score_;
	unsigned int DPMatrixDimensionRow;
	unsigned int DPMatrixDimensionColumn;
	// holds the aligned sequences
	std::string SeqA_align;
	std::string SeqB_align;
	std::string Symbol_align;
	std::unordered_map<int, int> Nucleotide_matching;
	// defines the DP matrices to be filled
	std::vector< std::vector<ALIGNSCORETYPE> > DPScoreMatrix;
	std::vector< std::vector<ALIGNSCORETYPE> > DPScoreMatrix_DEL;
	std::vector< std::vector<ALIGNSCORETYPE> > DPScoreMatrix_INS;
	std::vector< std::vector<enum AlignPath> > DPTraceMatrix;
	bool trace_back_called;
	int final_alignment_query_begin;
	int final_alignment_query_end;
	int final_alignment_target_begin;
	int final_alignment_target_end;
	// holds the scoring function
	ScoringFunction<ALIGNSCORETYPE>* ScoreScheme;
	// marks the alignment mode
	enum AlignMode Mode;
	bool Banded;
	unsigned int DownwardBand, RightwardBand;
	virtual inline void init_alignment(void);
	virtual void init_basic_parameters(
		  std::string &in_seqA, std::string &in_seqB, 
		  ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, 
		  enum AlignMode in_mode
	);
	virtual void align(void);
	virtual inline ALIGNSCORETYPE* get_cell_address(
	    const unsigned int i, const unsigned int j, const enum DPScoreMatrixName name
	);
	virtual inline char* get_character(std::string &seq, const unsigned int i);
	virtual inline std::vector<unsigned int> get_band_intersection(void);
	virtual inline void fill_initial_value(void);
	virtual inline void fill_DPMatrix_row(const unsigned int row);
	virtual inline void fill_cell_DPScoreMatrix_DEL(const unsigned int row, const unsigned int column);
	virtual inline void fill_cell_DPScoreMatrix_INS(const unsigned int row, const unsigned int column);
	virtual inline void fill_cell_DPScoreMatrix(const unsigned int row, const unsigned int column);
	virtual inline void fill_cell_DPTraceMatrix(const unsigned int row, const unsigned int column);
	virtual inline std::pair<unsigned int, unsigned int> compute_banded_region(const unsigned int row);
	virtual inline void check_band_setup(const int min_range);
	virtual inline ALIGNSCORETYPE get_global_alignment_score_best(void);
	virtual void set_max_score_index(int& idx_i, int& idx_j);
};

template<typename ALIGNSCORETYPE> SeqAlign<ALIGNSCORETYPE>::SeqAlign()	{
  best_score_ = -_MAX_VALUE; 
	return;
}


template<typename ALIGNSCORETYPE> 
SeqAlign<ALIGNSCORETYPE>::SeqAlign(
	  std::string &in_seqA, std::string &in_seqB, 
	  ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, enum AlignMode in_mode
)	{
	
	assert(in_seqA.length() > 0 && in_seqB.length() > 0);
	best_score_ = -_MAX_VALUE; 
	DPMatrixDimensionRow = in_seqA.length() + 1;
	DPMatrixDimensionColumn = in_seqB.length() + 1;
	init_basic_parameters(in_seqA, in_seqB, in_score_scheme, in_mode);
	Banded = false;
	return;
}

template<typename ALIGNSCORETYPE> SeqAlign<ALIGNSCORETYPE>::SeqAlign(
	std::string &in_seqA, std::string &in_seqB, ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, enum AlignMode in_mode, 
	unsigned int in_downward_band, unsigned int in_rightward_band, bool to_check_setup)	{
	
	assert(in_seqA.length() > 0 && in_seqB.length() > 0);
	assert(in_downward_band > 0 && in_rightward_band > 0);
	best_score_ = -_MAX_VALUE; 
	init_basic_parameters(in_seqA, in_seqB, in_score_scheme, in_mode);
	DPMatrixDimensionRow = in_seqA.length() + 1;
	DPMatrixDimensionColumn = in_seqB.length() + 1;
	Banded = true;
	DownwardBand = in_downward_band;
	RightwardBand = in_rightward_band;
	// for SEMIGLOBAL alignment, we have to guarantee seqA's full sequence is covered by the band
	if(Mode == SEMIGLOBAL && DownwardBand + SeqB.length() < SeqA.length())  {
	  DownwardBand = SeqA.length() - SeqB.length() + 1;
	}
	// checks the setup of the bands
	if(to_check_setup)	{
		int k = _min<int>(in_downward_band, in_rightward_band);
		check_band_setup(k);
	}
	return;
}

template<typename ALIGNSCORETYPE> SeqAlign<ALIGNSCORETYPE>::~SeqAlign()	{
	return;
}

template<typename ALIGNSCORETYPE>
void SeqAlign<ALIGNSCORETYPE>::SetMode(AlignMode in_mode)  {
  Mode = in_mode;
  //std::cout << "Set mode finished" << std::endl;
  return;
}

template<typename ALIGNSCORETYPE>
ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::GetBestScore() {
  return best_score_;
}

template<typename ALIGNSCORETYPE>
ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::GetBestGlobalScore()  {
  if(Mode == SEMIGLOBAL)  {
    // get the best score of the last row
    std::vector<unsigned int> intersection;
	  intersection = GetBandIntersection();
	  unsigned int i;
	  ALIGNSCORETYPE max_score = -_MAX_VALUE;
	  for(i = intersection[0]; i <= intersection[1]; ++ i)	{
		  if(*get_cell_address(DPMatrixDimensionRow - 1, i, FULL) > max_score)	{
		  	max_score = *get_cell_address(DPMatrixDimensionRow - 1, i, FULL);
		  }
	  }
	  return max_score;
  } else {
    return *get_cell_address(DPMatrixDimensionRow - 1, DPMatrixDimensionColumn - 1, FULL);
  }
}

template<typename ALIGNSCORETYPE>
std::vector<int> SeqAlign<ALIGNSCORETYPE>::GetAlignmentRegion(void){
  assert(trace_back_called);
  std::vector<int> region = 
      {final_alignment_query_begin, final_alignment_query_end, 
       final_alignment_target_begin, final_alignment_target_end};
  return region;
}

template<typename ALIGNSCORETYPE>
void SeqAlign<ALIGNSCORETYPE>::set_max_score_index(int& idx_i, int& idx_j)  {
  if(Mode == SEMILOCAL)  {
    unsigned int row_start = 0, column_start = 0;
    if(Banded)  {
      std::vector<unsigned int> intersection = get_band_intersection();
      row_start = intersection[0];
      column_start = intersection[2];
    }
    ALIGNSCORETYPE max_score = -_MAX_VALUE;
    unsigned int i;
    for(i = row_start; i < DPMatrixDimensionColumn; ++ i) {
      if(*get_cell_address(DPMatrixDimensionRow - 1, i, FULL) > max_score)  {
        max_score = *get_cell_address(DPMatrixDimensionRow - 1, i, FULL);
        idx_i = DPMatrixDimensionRow - 1;
        idx_j = i;
      }
    }
    for(i = column_start; i < DPMatrixDimensionRow; ++ i) {
      if(*get_cell_address(i, DPMatrixDimensionColumn - 1, FULL) > max_score)  {
        max_score = *get_cell_address(i, DPMatrixDimensionColumn - 1, FULL);
        idx_i = i;
        idx_j = DPMatrixDimensionColumn - 1;
      }
    }
    return;
  } else if(Mode == SEMIGLOBAL) {
    unsigned int row_start = 0;
    if(Banded)  {
      std::vector<unsigned int> intersection = get_band_intersection();
      row_start = intersection[0];
    }
    ALIGNSCORETYPE max_score = -_MAX_VALUE;
    unsigned int i;
    for(i = row_start; i < DPMatrixDimensionColumn; ++ i) {
      if(*get_cell_address(DPMatrixDimensionRow - 1, i, FULL) > max_score)  {
        max_score = *get_cell_address(DPMatrixDimensionRow - 1, i, FULL);
        idx_i = DPMatrixDimensionRow - 1;
        idx_j = i;
      }
    }
    //std::cout << idx_i << " " << idx_j << std::endl;
    return;
  } else if(Mode == GLOBAL) {
    idx_i = DPMatrixDimensionRow - 1;
    idx_j = DPMatrixDimensionColumn - 1;
    return;
  } else if(Mode == LOCAL) {
    ALIGNSCORETYPE max_score = -_MAX_VALUE;
    unsigned int i, j;
    for(i = 0; i < DPMatrixDimensionRow; ++ i) {
	    unsigned int left_bound = 1, right_bound = DPMatrixDimensionColumn - 1;
	    if(Banded)	{
		    std::pair<unsigned int, unsigned int> bound = compute_banded_region(i);
		    left_bound = bound.first;
		    right_bound = bound.second > DPMatrixDimensionColumn - 1 ? DPMatrixDimensionColumn - 1 : bound.second;
	    }
	    for(j = left_bound; j <= right_bound; ++ j)	{
        if(*get_cell_address(i, j, FULL) > max_score) {
          max_score = *get_cell_address(i, j, FULL);
          idx_i = i;
          idx_j = j;
        }
	    }
	  }
	  return;
  }
}

template<typename ALIGNSCORETYPE> 
ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::ComputeHammingDistance(void) {
  assert(SeqA.length() == SeqB.length());
  ALIGNSCORETYPE score = 0;
  for(unsigned int i = 0; i < SeqA.length(); ++ i) {
    score += ScoreScheme->CheckMatchScore(SeqA[i], SeqB[i]);
  }
  return score;
}

/* TODO: get a more sophisticated GetAlignmentScore function depending on the banded intersection region */
template<typename ALIGNSCORETYPE>
inline ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::GetAlignmentScore()	{
	ALIGNSCORETYPE optimal_score = -_MAX_VALUE;
	if(Mode == GLOBAL)	{
		optimal_score = get_global_alignment_score_best();
	}
	return optimal_score;
}

template<typename ALIGNSCORETYPE>
void SeqAlign<ALIGNSCORETYPE>::FillEdgeScore(EdgeScoreBundle<ALIGNSCORETYPE>& edge_score) {
  edge_score.row = GetRowScore(DPMatrixDimensionRow - 1, FULL);
  edge_score.row_del = GetRowScore(DPMatrixDimensionRow - 1, DEL);
  edge_score.row_ins = GetRowScore(DPMatrixDimensionRow - 1, INS);
  edge_score.column = GetColumnScore(DPMatrixDimensionColumn - 1, FULL);
  edge_score.column_del = GetColumnScore(DPMatrixDimensionColumn - 1, DEL);
  edge_score.column_ins = GetColumnScore(DPMatrixDimensionColumn - 1, INS);
  std::vector<unsigned int> intersection = GetBandIntersection();
  edge_score.begin_row = intersection[0];
  edge_score.begin_column = intersection[2];
  return;
}

template<typename ALIGNSCORETYPE>
inline std::vector<ALIGNSCORETYPE> SeqAlign<ALIGNSCORETYPE>::GetRowScore(unsigned int row, enum DPScoreMatrixName name)	{
	assert(row >= 0 && row < DPMatrixDimensionRow);
	std::vector<ALIGNSCORETYPE> row_score(DPMatrixDimensionColumn);
	//std::cout << "GetRowScore called!!!" << std::endl;
	
	for(unsigned int i = 0; i < DPMatrixDimensionColumn; ++ i)	{
		if(name == TRACE)	{
			row_score[i] = *(enum AlignPath*)get_cell_address(row, i, name);
		}	else	{
			row_score[i] = *get_cell_address(row, i, name);
		}
		//std::cout << "-	" << row_score[i] << std::endl;
	}
	return row_score;
}

template<typename ALIGNSCORETYPE>
inline std::vector<ALIGNSCORETYPE> SeqAlign<ALIGNSCORETYPE>::GetColumnScore(unsigned int column, enum DPScoreMatrixName name)	{
	assert(column >= 0 && column < DPMatrixDimensionColumn);
	std::vector<ALIGNSCORETYPE> column_score(DPMatrixDimensionRow);
	//std::cout << "GetColumnScore called!!!" << std::endl;
	for(unsigned int i = 0; i < DPMatrixDimensionRow; ++ i)	{
		if(name == TRACE)	{
			column_score[i] = *(enum AlignPath*)get_cell_address(i, column, name);
		}	else	{
			column_score[i] = *get_cell_address(i, column, name);
		}
		//std::cout << "|	" << column_score[i] << std::endl;
	}
	return column_score;
}

template<typename ALIGNSCORETYPE>
inline std::vector<unsigned int> SeqAlign<ALIGNSCORETYPE>::GetBandIntersection()	{
	return get_band_intersection();
}

template<typename ALIGNSCORETYPE>
inline ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::GetBestEdgeScore()	{
	std::vector<unsigned int> intersection;
	intersection = GetBandIntersection();
	unsigned int i;
	ALIGNSCORETYPE max_score = -_MAX_VALUE;
	for(i = intersection[0]; i <= intersection[1]; ++ i)	{
		if(*get_cell_address(DPMatrixDimensionRow - 1, i, FULL) > max_score)	{
			max_score = *get_cell_address(DPMatrixDimensionRow - 1, i, FULL);
		}
	}
	for(i = intersection[2]; i <= intersection[3]; ++ i)	{
		if(*get_cell_address(i, DPMatrixDimensionColumn - 1, FULL) > max_score)	{
			max_score = *get_cell_address(i, DPMatrixDimensionColumn - 1, FULL);
		}
	}
	return max_score;
}

template<typename ALIGNSCORETYPE>
inline ALIGNSCORETYPE SeqAlign<ALIGNSCORETYPE>::get_global_alignment_score_best()	{
	assert(Mode == GLOBAL);
	if(Banded)	{
		std::vector<unsigned int> intersection(4);
		intersection = get_band_intersection();
		assert(intersection[0] < DPMatrixDimensionColumn);
		assert(intersection[2] < DPMatrixDimensionRow);
	}
	return *get_cell_address(DPMatrixDimensionRow - 1, DPMatrixDimensionColumn - 1, FULL);
}

template<typename ALIGNSCORETYPE>
inline std::vector<unsigned int> SeqAlign<ALIGNSCORETYPE>::get_band_intersection(void)	{
	// computing the intersection between the band and the two axies
	// intersection[0] and intersection[1]: intersection coordinates in the x-axis
	// intersection[2] and intersection[3]: intersection coordinates in the y-axis
	assert(Banded);
	assert(DownwardBand > 0 && RightwardBand > 0);
	std::vector<unsigned int> intersection(4);
	//std::cout << "In get_band_intersection:  " << DPMatrixDimensionRow << "  " << DPMatrixDimensionColumn << " " << DownwardBand << "  " << RightwardBand << std::endl;
	intersection[0] = DPMatrixDimensionRow - 1 - DownwardBand + 1;
	intersection[1] = DPMatrixDimensionRow - 1 + RightwardBand - 1;
	intersection[2] = DPMatrixDimensionColumn - 1 - RightwardBand + 1;
	intersection[3] = DPMatrixDimensionColumn - 1 + DownwardBand - 1;
	
	intersection[0] = static_cast<int>(intersection[0]) < 0 ? 0 : intersection[0];
	intersection[1] = intersection[1] > DPMatrixDimensionColumn - 1 ? DPMatrixDimensionColumn - 1 : intersection[1];
	intersection[2] = static_cast<int>(intersection[2]) < 0 ? 0 : intersection[2];
	intersection[3] = intersection[3] > DPMatrixDimensionRow - 1 ? DPMatrixDimensionRow - 1 : intersection[3];
	
	return intersection;
}

template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::init_basic_parameters(
	std::string &in_seqA, std::string &in_seqB, ScoringFunction<ALIGNSCORETYPE>* in_score_scheme, enum AlignMode in_mode)	{
	assert(in_seqA.length() > 0 && in_seqB.length() > 0);
	SeqA = in_seqA;
	SeqB = in_seqB;
	ScoreScheme = in_score_scheme;
	Mode = in_mode;
	return;
}


template<typename ALIGNSCORETYPE> 
inline char* SeqAlign<ALIGNSCORETYPE>::get_character(std::string &seq, const unsigned int i)	{
	assert(i >= 0 && i < seq.length());
	return &seq[i];
}


template<typename ALIGNSCORETYPE> ALIGNSCORETYPE* SeqAlign<ALIGNSCORETYPE>::get_cell_address(const unsigned int i, const unsigned int j, const enum DPScoreMatrixName name)	{
  //std::cout << "problem of get_cell_address: " << i << " " <<  j << std::endl;
  if(!(i >= 0 && j >= 0 && i < DPMatrixDimensionRow && j < DPMatrixDimensionColumn))  {
    //std::cout << SeqA << "  " << SeqB << std::endl;
    std::cout << "problem of get_cell_address: " << i << " " <<  j << std::endl;
    exit(0);
  }
	assert(i >= 0 && j >= 0 && i < DPMatrixDimensionRow && j < DPMatrixDimensionColumn);
	if(name == FULL)	{
		return &DPScoreMatrix[i][j];
	}	else if(name == DEL)	{
		return &DPScoreMatrix_DEL[i][j];
	}	else if(name == INS)	{
		return &DPScoreMatrix_INS[i][j];
	}	else if(name == TRACE)	{
		return (ALIGNSCORETYPE*) &DPTraceMatrix[i][j];
	}
	return NULL;
}


// Computes the banded region for a given row
template<typename ALIGNSCORETYPE> 
std::pair<unsigned int, unsigned int> SeqAlign<ALIGNSCORETYPE>::compute_banded_region(const unsigned int row)	{
	assert(Banded);
	assert(DownwardBand > 0);
	assert(RightwardBand > 0);
	assert(row >= 0 && row < DPMatrixDimensionRow);
	// computes the corresponding banded region
	std::pair<unsigned int, unsigned int> region;
	if(row <= DownwardBand)	{
		region.first = 1;
	}	else	{
		region.first = row - DownwardBand + 1;
	}	
	region.second = row + RightwardBand - 1;
	return region;
}

// Checks the banding setup
template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::check_band_setup(const int min_range)	{
	assert(Banded);
	// preset the minimum range that should be allowed for the banded alignment
	if(DownwardBand >= DPMatrixDimensionRow && RightwardBand >= DPMatrixDimensionColumn)	{
		std::cout << "Warning: band setup is too loose. Banded alignment mode disabled." << std::endl;
		Banded = false;
		return;
	}	else if(DownwardBand >= DPMatrixDimensionRow)	{
		std::cout << "Warning: downward band setup is too loose. Auto adjust downward band to sequence length." << std::endl;
		DownwardBand = DPMatrixDimensionRow;
	}	else if(RightwardBand >= DPMatrixDimensionColumn)	{
		std::cout << "Warning: rightward band setup is too loose. Auto adjust rightward band to sequence length." << std::endl;
		RightwardBand = DPMatrixDimensionColumn;
	}
	// compute corresponding downward and rightward band setup
	unsigned int min_downward_band = min_range + SeqA.length() - SeqB.length();
	unsigned int min_rightward_band = min_range + SeqB.length() - SeqA.length();
	if(DownwardBand < min_downward_band)	{
		DownwardBand = min_downward_band;
		std::cout << "Warning: band setup is too strict. Auto adjust downward band to " << DownwardBand << "." << std::endl;
	}
	if(RightwardBand < min_rightward_band)	{
		RightwardBand = min_rightward_band;
		std::cout << "Warning: band setup is too strict. Auto adjust rightward band to " << RightwardBand << "." << std::endl;
	}
	return;
}

// Initilize DP matrix space for alignment
template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::init_alignment(void)	{
	
	unsigned int m = DPMatrixDimensionRow;
	unsigned int n = DPMatrixDimensionColumn;
	assert(m > 1 && n > 1);
	
	DPScoreMatrix.resize(m);
	DPScoreMatrix_INS.resize(m);
	DPScoreMatrix_DEL.resize(m);
	DPTraceMatrix.resize(m);
	for(unsigned i = 0; i < m; ++ i)	{
		DPScoreMatrix[i].resize(n, -_MAX_VALUE);
		DPScoreMatrix_INS[i].resize(n, -_MAX_VALUE);
		DPScoreMatrix_DEL[i].resize(n, -_MAX_VALUE);
		DPTraceMatrix[i].resize(n);
	}
	
	return;
}

// Parsing arguments and prepares alignment
template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::Align(void)	{
	//std::cout << "before init alignment" << std::endl;
	//if(Mode == SEMIGLOBAL)  {
	  //std::cout << SeqA << "  " << SeqB << std::endl;
	//}
	init_alignment();
	//std::cout << "after init alignment" << std::endl;
	align();
	return;
}

// Fills the boundary values
template<typename ALIGNSCORETYPE> inline void SeqAlign<ALIGNSCORETYPE>::fill_initial_value(void)	{
	assert(SeqA.length() > 0 && SeqB.length() > 0);
	unsigned int i;
	*get_cell_address(0, 0, FULL) = ScoreScheme->GetGapOpen();
	
	// determine the boundary to be filled
	unsigned int down_bound = 
	    (Banded && DownwardBand - 1 < DPMatrixDimensionRow) ? DownwardBand - 1 : DPMatrixDimensionRow - 1;
	unsigned int right_bound = 
	    (Banded && RightwardBand - 1 < DPMatrixDimensionColumn) ? RightwardBand - 1 : DPMatrixDimensionColumn - 1;
	//std::cout << "fill_initial_value:: dimensions: " << DPMatrixDimensionRow << "  " << DPMatrixDimensionColumn << std::endl; 
	//std::cout << "fill_initial_value:: down and right bound: " << down_bound << "  " << right_bound << std::endl; 
	for(i = 1; i <= down_bound; ++ i)	{
	  if(Mode == GLOBAL || Mode == SEMIGLOBAL)  {
		  *get_cell_address(i, 0, FULL) = *get_cell_address(i - 1, 0, FULL) + ScoreScheme->GetGapExtend();
		  *get_cell_address(i, 0, DEL) = *get_cell_address(i, 0, FULL) + ScoreScheme->GetGapOpen();
		  *get_cell_address(i, 0, INS) = *get_cell_address(i, 0, FULL);
	  	*(enum AlignPath*)get_cell_address(i, 0, TRACE) = UP;
		} else if(Mode == LOCAL || Mode == SEMILOCAL)  {
		  *get_cell_address(i, 0, FULL) = *get_cell_address(i, 0, DEL) = *get_cell_address(i, 0, INS) = 0;
		  *(enum AlignPath*)get_cell_address(i, 0, TRACE) = TERMINAL;
		}
	}
	for(i = 1; i <= right_bound; ++ i)	{
	  if(Mode == GLOBAL)  {
		  *get_cell_address(0, i, FULL) = *get_cell_address(0, i - 1, FULL) + ScoreScheme->GetGapExtend();
		  *get_cell_address(0, i, INS) = *get_cell_address(0, i, FULL) + ScoreScheme->GetGapOpen();
	  	*get_cell_address(0, i, DEL) = *get_cell_address(0, i, FULL);
	  	// if the mode is semiglobal then we do not need to go left any more
		  *(enum AlignPath*)get_cell_address(0, i, TRACE) = LEFT;
		} else if(Mode == LOCAL || Mode == SEMILOCAL || Mode == SEMIGLOBAL) {
		  *get_cell_address(0, i, FULL) = *get_cell_address(0, i, DEL) = *get_cell_address(0, i, INS) = 0;
		  *(enum AlignPath*)get_cell_address(0, i, TRACE) = TERMINAL;
		}
	}
	*get_cell_address(0, 0, FULL) = *get_cell_address(0, 0, DEL) = *get_cell_address(0, 0, INS) = 0;
	*get_cell_address(0, 0, TRACE) = TERMINAL;
	//std::cout << "fill_initial_value:: good after filling boundaries" << std::endl;
	if(Banded)	{
		// filling special values
		if(DownwardBand < DPMatrixDimensionRow)	{
			*get_cell_address(DownwardBand, 0, FULL) = *get_cell_address(DownwardBand, 0, DEL) = -_MAX_VALUE;
		}
		if(RightwardBand < DPMatrixDimensionColumn)	{
			*get_cell_address(0, RightwardBand, FULL) = *get_cell_address(0, RightwardBand, INS) = -_MAX_VALUE;
		}
		for(i = 1; i < DPMatrixDimensionRow; ++ i)	{
			std::pair <unsigned int, unsigned int> bound = compute_banded_region(i);
			if(bound.first - 1 > 0 && bound.first - 1 < DPMatrixDimensionColumn)	{
				*get_cell_address(i, bound.first - 1, FULL) = *get_cell_address(i, bound.first - 1, DEL) = -_MAX_VALUE;
			}
			if(bound.second + 1 > 0 && bound.second + 1 < DPMatrixDimensionColumn)	{
				*get_cell_address(i, bound.second + 1, FULL) = *get_cell_address(i, bound.second + 1, INS) = -_MAX_VALUE;
			}
		}
	}
	//std::cout << "fill_initial_value:: good before return" << std::endl;
	return;
}

// Fills a cell in the DPScoreMatrix
template<typename ALIGNSCORETYPE> inline void SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix(const unsigned int row, const unsigned int column)	{
	assert(row > 0 && row < DPMatrixDimensionRow && column > 0 && column < DPMatrixDimensionColumn);
	
	*get_cell_address(row, column, FULL) = _max<ALIGNSCORETYPE>(
		*get_cell_address(row, column, DEL), *get_cell_address(row, column, INS)
	);
	*get_cell_address(row, column, FULL) = _max<ALIGNSCORETYPE>(
		*get_cell_address(row, column, FULL),
		*get_cell_address(row - 1, column - 1, FULL) + 
			ScoreScheme->CheckMatchScore(*get_character(SeqA, row - 1), *get_character(SeqB, column - 1))
	);
	if(Mode == LOCAL && *get_cell_address(row, column, FULL) <= 0)  {
	  *get_cell_address(row, column, FULL) = 0;
	}	  	
	if(*get_cell_address(row, column, FULL) > best_score_)  {
	  best_score_ = *get_cell_address(row, column, FULL);
	}
	return;
}

// Fills a cell in the DPTraceMatrix
template<typename ALIGNSCORETYPE> 
inline void SeqAlign<ALIGNSCORETYPE>::fill_cell_DPTraceMatrix(const unsigned int row, const unsigned int column)	{
	assert(row > 0 && row < DPMatrixDimensionRow && column > 0 && column < DPMatrixDimensionColumn);
	if(Mode == LOCAL && *get_cell_address(row, column, FULL) <= 0)  {
	  *(enum AlignPath*)get_cell_address(row, column, TRACE) = TERMINAL;
	} else if(*get_cell_address(row, column, FULL) <= *get_cell_address(row, column, DEL))	  {
		*(enum AlignPath*)get_cell_address(row, column, TRACE) = LEFT;
	}	else if(*get_cell_address(row, column, FULL) <= *get_cell_address(row, column, INS))	  {
		*(enum AlignPath*)get_cell_address(row, column, TRACE) = UP;
	}	else	{
		*(enum AlignPath*)get_cell_address(row, column, TRACE) = DIAGONAL;
	}
	return;
}

// Fillls a cell in the DPScoreMatrix_DEL
template<typename ALIGNSCORETYPE> 
inline void SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix_DEL(const unsigned int row, const unsigned int column)	{
	assert(row > 0 && row < DPMatrixDimensionRow && column > 0 && column < DPMatrixDimensionColumn);
	*get_cell_address(row, column, DEL) = _max<ALIGNSCORETYPE>(
		*get_cell_address(row, column - 1, DEL) + ScoreScheme->GetGapExtend(),
		*get_cell_address(row, column - 1, FULL) + ScoreScheme->GetGapOpen() + ScoreScheme->GetGapExtend()
	);
	if(Mode == LOCAL && *get_cell_address(row, column, DEL) <= 0)  {
	  *get_cell_address(row, column, DEL) = 0;
	}
	return;
}

// Fills a cell in the DPScoreMatrix_INS
template<typename ALIGNSCORETYPE> 
inline void SeqAlign<ALIGNSCORETYPE>::fill_cell_DPScoreMatrix_INS(const unsigned int row, const unsigned int column)	{
	assert(row > 0 && row < DPMatrixDimensionRow && column > 0 && column < DPMatrixDimensionColumn);
	*get_cell_address(row, column, INS) = _max<ALIGNSCORETYPE>(
		*get_cell_address(row - 1, column, INS) + ScoreScheme->GetGapExtend(),
		*get_cell_address(row - 1, column, FULL) + ScoreScheme->GetGapOpen() + ScoreScheme->GetGapExtend()
	);
	if(Mode == LOCAL && *get_cell_address(row, column, INS) <= 0)  {
	  *get_cell_address(row, column, INS) = 0;
	}
	return;
}


// Fills a row of the DP matrices
template<typename ALIGNSCORETYPE> 
inline void SeqAlign<ALIGNSCORETYPE>::fill_DPMatrix_row(const unsigned int row)	{
	assert(row > 0 && row < DPMatrixDimensionRow);
	// check if banded alignment is needed
	unsigned int left_bound = 1, right_bound = DPMatrixDimensionColumn - 1;
	if(Banded)	{
		std::pair<unsigned int, unsigned int> bound = compute_banded_region(row);
		left_bound = bound.first;
		right_bound = bound.second > DPMatrixDimensionColumn - 1 ? DPMatrixDimensionColumn - 1 : bound.second;
	}
	unsigned int i;
	for(i = left_bound; i <= right_bound; ++ i)	{
		fill_cell_DPScoreMatrix_DEL(row, i);
		fill_cell_DPScoreMatrix_INS(row, i);
		fill_cell_DPScoreMatrix(row, i);
		fill_cell_DPTraceMatrix(row, i);
	}
	return;
}	



// Trace back function to retrive the alignment
template<typename ALIGNSCORETYPE> 
void SeqAlign<ALIGNSCORETYPE>::TraceBack(void)	{
	trace_back_called = true;
	if(Banded)	{
		std::vector<unsigned int> intersection(4);
		intersection = get_band_intersection();
		//assert(intersection[0] < DPMatrixDimensionColumn);
		//assert(intersection[2] < DPMatrixDimensionRow);
	}
	int idx_i, idx_j;
	set_max_score_index(idx_i, idx_j);
	final_alignment_query_end = final_alignment_query_begin = idx_i - 1;
	final_alignment_target_end = final_alignment_target_begin = idx_j - 1;
	alignment_score = *get_cell_address(idx_i, idx_j, FULL);
	//std::cout << "TraceBack:: indices: " << idx_i << " " << idx_j << std::endl;
	while((idx_i > 0 || idx_j > 0) && *(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) != TERMINAL)	{
	  if(Mode == SEMIGLOBAL && idx_i == 0)  {
	    break;
	  }
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


// Core algorithm that does the alignment
template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::align(void)	{
	
	assert(SeqA.length() > 0 && SeqB.length() > 0);	
	fill_initial_value();
	//std::cout << "align::finish filling initial values" << std::endl;
	unsigned int i;	
	for(i = 1; i < DPMatrixDimensionRow; ++ i)	{
		fill_DPMatrix_row(i);
	}
	/* TODO: write a global wrapper function TraceBack, which determines which traceback
		route to call */
	//trace_back();
	return;
}	

template<typename ALIGNSCORETYPE>
void SeqAlign<ALIGNSCORETYPE>::ListAllHSPPositions(
    ALIGNSCORETYPE score_cutoff, 
    std::list<std::vector<int> >& hsp_range
)  {
  // ensure that the current setting is local alignment
  assert(Mode == LOCAL);
  // setting an indicator array to check if the cell has been visited before
  int i, j;
  bool **visit_indicator = new bool*[DPMatrixDimensionRow];
  for(i = 0; i < (int) DPMatrixDimensionRow; ++ i) {
    visit_indicator[i] = new bool[DPMatrixDimensionColumn];
    for(j = 0; j < (int) DPMatrixDimensionColumn; ++ j) {
      visit_indicator[i][j] = false;
    }
  }

  // finds all alignments
  for(i = DPMatrixDimensionRow - 1; i >= 0; -- i) {
    for(j = DPMatrixDimensionColumn - 1; j >= 0; -- j) {
      // check if the cell has been visited
      //std::cout << "index outside while loop: " << i << " " << j << std::endl;
      if(visit_indicator[i][j])  {
        continue;
      }
      // otherwise if the cell is not visited
      if(*get_cell_address(i, j, FULL) >= score_cutoff)  {
        int max_i = i, max_j = j;
        int idx_i = i, idx_j = j;
        bool has_conflict = false;  // indicate whether the path has some overlap with existing paths
        std::string hsp_alignment_A, hsp_alignment_B;
        while((idx_i > 0 || idx_j > 0) && *(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) != TERMINAL)	{
          //std::cout << "index inside while: " << idx_i << " " << idx_j << std::endl;
          if(*get_cell_address(idx_i, idx_j, FULL) > *get_cell_address(max_i, max_j, FULL))  {
            // setting the maximum entry
            max_i = idx_i;
            max_j = idx_j;
            // dropping the existing alignment
            hsp_alignment_A = "";
            hsp_alignment_B = "";
          }
          if(!visit_indicator[idx_i][idx_j])  {
            visit_indicator[idx_i][idx_j] = true;
          } else  {
            has_conflict = true;
            break;
          }
      		if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == DIAGONAL)	{
			      -- idx_i;
			      -- idx_j;
			      hsp_alignment_A += SeqA[idx_i];
			      hsp_alignment_B += SeqB[idx_j];
      		}	else if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == LEFT)	{
			      -- idx_j;
      			hsp_alignment_A += '-';
			      hsp_alignment_B += SeqB[idx_j];
      		}	else if(*(enum AlignPath*)get_cell_address(idx_i, idx_j, TRACE) == UP)	{
			      -- idx_i;
			      hsp_alignment_A += SeqA[idx_i];
			      hsp_alignment_B += '-';
      		}
	      } // end while
	      if(!has_conflict || *get_cell_address(i, j, FULL) - *get_cell_address(idx_i, idx_j, FULL) > score_cutoff)  {
	        // record the alignment
	        std::vector<int> hsp_coordinates(4);
	        hsp_coordinates = {max_i, max_j, idx_i, idx_j};
	        hsp_range.push_back(hsp_coordinates);
	      }
      } // end if
    }
  }
  // collect memory
  for(i = 0; i < (int) DPMatrixDimensionRow; ++ i) {
    delete [] visit_indicator[i];
  }
  delete [] visit_indicator;
  return;
}

// Printing the alignment
template<typename ALIGNSCORETYPE> void SeqAlign<ALIGNSCORETYPE>::PrintAlignment(void)	{
	unsigned int i;
	double bit_score = ScoreScheme->ComputeBitScore(alignment_score);
	std::cout << "bit score:  " << bit_score << "(" << alignment_score << ")" << std::endl;
	for(i = 0; i < SeqA_align.length(); i += 100)	{
		std::cout << "\t" << SeqA_align.substr(i, 100) << std::endl;
		std::cout << "\t" << Symbol_align.substr(i, 100) << std::endl;
		std::cout << "\t" << SeqB_align.substr(i, 100) << std::endl;
		std::cout << std::endl;
	}
	return;	
}

template<typename ALIGNSCORETYPE> 
void SeqAlign<ALIGNSCORETYPE>::GetAlignment(
    std::unordered_map<int, int>& nuc_match, 
    std::string& seqA_align, std::string& seqB_align, std::string& symbol_align
) {
  nuc_match = Nucleotide_matching;
  seqA_align = SeqA_align;
  seqB_align = SeqB_align;
  symbol_align = Symbol_align;
  return;
}

#endif
