#include "hmm_profile.h"
#include "seq_align.h"

#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <tuple>

#ifndef _SEQ_ALIGN_EXTEND_P_
#define _SEQ_ALIGN_EXTEND_P_

typedef int SCORETYPE;
enum DIRECTION {FORWARD = 0, BACKWARD = 1};
struct ScoreRecordType  {
  SCORETYPE *mmr, *mmc; // match matrix row and column
  SCORETYPE *imr, *imc; // insert matrix row and column
  SCORETYPE *dmr, *dmc; // delete matrix row and column
  int r_size, c_size; // sizes of row and columns
  int rbegin, cbegin; // begins of the scores (accoding to band settings)
};
struct AlnResultStruct  {
  // note that m_score - n_score corresponds to the minimum (optimal) score can be achieved
  SCORETYPE m_score;    // score computed from the model
  SCORETYPE n_score;    // score computed from null hypothesis
  int m_index;          // index in the model that achieves the optima
  int s_index;          // index in the sequence that achieves the optima
  int m_filled;         // the row-index where alignment is computed
  int s_filled;         // the column-index where alignment is computed
};

// Extending profile-sequence alignment (using banded Viterbi) assuming an anchor is found.
// Tries to extend to the end of the profile.
class SeqAlignExtendP {
 public:
  SeqAlignExtendP();
  SeqAlignExtendP(HMMProfile &in_model, int in_band_size);
  ~SeqAlignExtendP();
  void AlignInitFW(int mbegin, std::string &seq, AlnResultStruct &result);
  void AlignExtdFW(
    int mbegin,    
    int m_filled, int s_filled,
    std::string &seq,
    AlnResultStruct &result
  );
  void AlignInitRE(int mbegin, std::string &seq, AlnResultStruct &result);
  void AlignExtdRE(
    int mbegin,    
    int m_filled, int s_filled,
    std::string &seq,
    AlnResultStruct &result
  );
  friend class AssembleExtendP;
 protected: 
  bool is_model_set_, is_band_set_, is_matrix_initialized_;
  HMMProfile model_;
  int pflen_;       // the length of the HMM profile, fixed for each instance
  int band_size_;   // band size for DP
  SCORETYPE **mmx_; // scoring matrix assuming the last operation is a match   
  SCORETYPE **dmx_; // scoring matrix assuming the last operation is a deletion
  SCORETYPE **imx_; // scoring matrix assuming the last operation is an insertion
  // null hypothesis parsing score for the sequence
  // use vector for resizing simplicity
  std::vector<SCORETYPE> null_score_;    
  
  // compute the beginning and ending index of the valid region in the row
  inline std::pair<int, int> GetBandRegionRow(int r);
  // compute the score for cell [i][j] in the matrices, forward direction
  // mbegin: head index in HMM model whose suffix is to be aligned
  // seq_shift: length of sequence that have already been aligned for the extension (0 for init)
  inline void FillScoreCoreFW(int i, int j, int mbegin, std::string &seq, int seq_shift);
  // compute the score for cell [i][j] in the matrices, reverse direction
  // mbegin: head index in HMM model whose suffix is to be aligned
  // seq_shift: length of sequence that have already been aligned for the extension (0 for init)
  inline void FillScoreCoreRE(int i, int j, int mbegin, std::string &seq, int seq_shift);
  void AllocMatrix(void);
  void InitMatrix(void);
  void CalNull(int pivot, std::string &seq);
  
};

inline int Min(int a, int b)  {
  return a >= b ? b : a;
}

inline int Max(int a, int b)  {
  return a >= b ? a : b;
}

inline std::pair<int, int> SeqAlignExtendP::GetBandRegionRow(int r)	{
	// computes the corresponding banded region
	std::pair<int, int> region;
	if(r <= band_size_)	{
		region.first = 0;
	}	else	{
		region.first = r - band_size_;
	}	
	region.second = r + band_size_;
	return region;
}

inline void SeqAlignExtendP::FillScoreCoreFW(
    int i, int j, 
    int mbegin, std::string &seq, 
    int seq_shift
)  {
  mmx_[i][j] = Min(
      imx_[i - 1][j - 1] + model_.TranScore(i + mbegin - 2, IM), 
      dmx_[i - 1][j - 1] + model_.TranScore(i + mbegin - 2, DM)
  );
  mmx_[i][j] = Min(
      mmx_[i - 1][j - 1] + model_.TranScore(i + mbegin - 2, MM),
      mmx_[i][j]
  );
  // include the emission probability
  mmx_[i][j] += model_.MEmitScore(i + mbegin - 1, seq[j - seq_shift - 1]);  
  // Compute insert state
  imx_[i][j] = model_.IEmitScore(i + mbegin - 1, seq[j - seq_shift - 1]) + Min(
      mmx_[i][j - 1] + model_.TranScore(i + mbegin - 1, MI),
      imx_[i][j - 1] + model_.TranScore(i + mbegin - 1, II)
  );
  // Compute delete state, note that there is no emission for delete states
  dmx_[i][j] = Min(
      mmx_[i - 1][j] + model_.TranScore(i + mbegin - 2, MD),
      dmx_[i - 1][j] + model_.TranScore(i + mbegin - 2, DD)
  );
  return;
}

inline void SeqAlignExtendP::FillScoreCoreRE(
    int i, int j, 
    int mbegin, std::string &seq, 
    int seq_shift
)  {
  mmx_[i][j] = Min(
      imx_[i - 1][j - 1] + model_.TranScore(mbegin - i + 1, MI), 
      dmx_[i - 1][j - 1] + model_.TranScore(mbegin - i + 1, MD)
  );
  mmx_[i][j] = Min(
      mmx_[i - 1][j - 1] + model_.TranScore(mbegin - i + 1, MM),
      mmx_[i][j]
  );
  //cout << mmx_[i][j] << " " << mmx_[i - 1][j - 1] << "  " << model_.TranScore(mbegin - i + 1, MM) << endl;
  // include the emission probability
  mmx_[i][j] += model_.MEmitScore(mbegin - i + 1, seq[j - seq_shift - 1]);  
  // Compute insert state
  imx_[i][j] = model_.IEmitScore(mbegin - i + 1, seq[j -  seq_shift - 1]) + Min(
      mmx_[i][j - 1] + model_.TranScore(mbegin - i + 1, IM),
      imx_[i][j - 1] + model_.TranScore(mbegin - i + 1, II)
  );
  // Compute delete state, note that there is no emission for delete states
  dmx_[i][j] = Min(
      mmx_[i - 1][j] + model_.TranScore(mbegin - i + 1, DM),
      dmx_[i - 1][j] + model_.TranScore(mbegin - i + 1, DD)
  );
  return;
}

#endif
