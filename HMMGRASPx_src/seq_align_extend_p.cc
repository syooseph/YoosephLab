#include "seq_align_extend_p.h"

using namespace std;


SeqAlignExtendP::SeqAlignExtendP() {
  is_model_set_ = is_band_set_ = is_matrix_initialized_ = false;
  return;
}

SeqAlignExtendP::SeqAlignExtendP(HMMProfile &in_model, int in_band_size)  {
  model_ = in_model; pflen_ = in_model.GetProfileLen();
  band_size_ = in_band_size;
  is_model_set_ = true;
  is_band_set_ = true;
  AllocMatrix();    // allocate matrix space accoding to band size and model length
  InitMatrix();     // initialize the matrices for boundary conditions
  // do not initialize null_score_ matrix until target sequence is present
  return;
}

SeqAlignExtendP::~SeqAlignExtendP() {
  if(is_matrix_initialized_)  {
    for(int i = 0; i < pflen_ + 1; ++ i) {
      delete [] mmx_[i];
      delete [] imx_[i];
      delete [] dmx_[i];
    }
    delete [] mmx_;
    delete [] imx_;
    delete [] dmx_;
  }
  return;
}

void SeqAlignExtendP::AllocMatrix(void) {
  // initialize the matrices
  mmx_ = new int* [pflen_ + 1];
  imx_ = new int* [pflen_ + 1];
  dmx_ = new int* [pflen_ + 1];
  for(int i = 0; i < pflen_ + 1; ++ i) {
    mmx_[i] = new int [pflen_ + band_size_ + 1];
    imx_[i] = new int [pflen_ + band_size_ + 1];
    dmx_[i] = new int [pflen_ + band_size_ + 1];
  }
  is_matrix_initialized_ = true;
  return;
}

void SeqAlignExtendP::InitMatrix(void)  {
  int i;
  // setting default value for each row
  int *dv = new int [pflen_ + band_size_ + 1];
  for(i = 0; i < pflen_ + band_size_ + 1; ++ i) {
    // note that we need to minimize the scores, which are negative log probabilities
    // in this case set the impossibles to positive infinity
    dv[i] = INF;  
  }
  for(i = 0; i < pflen_ + 1; ++ i) {
    memcpy(mmx_[i], dv, 4 * (pflen_ + band_size_ + 1));
    memcpy(imx_[i], dv, 4 * (pflen_ + band_size_ + 1));
    memcpy(dmx_[i], dv, 4 * (pflen_ + band_size_ + 1));
  }
  delete [] dv;
  return;
}

// add the null emission score the null_score_
// first resize the matrix to consider new sequence, then compute the score
// pivot is the place in the vector to insert null scor for seq
void SeqAlignExtendP::CalNull(int pivot, std::string &seq) {
  assert((unsigned int) pivot <= null_score_.size());
  unsigned int i;
  if(pivot == 0)  {
    // seq is the very first sequence to be considered
    null_score_.resize(seq.length());
    null_score_[0] = model_.MEmitScore(0, seq[0]);
    for(i = 1; i < seq.length(); ++ i) {
      null_score_[i] = null_score_[i - 1] + model_.MEmitScore(0, seq[i]);
    }
  } else  {
    // seq is an appended sequence
    null_score_.resize(pivot + seq.length());
    for(i = 0; i < seq.length(); ++ i) {
      null_score_[pivot + i] = null_score_[pivot + i - 1] + model_.MEmitScore(0, seq[i]);
    }
  }
  return;
}

// the very first alignment with a sequence intput, the entire
// sequence is to be aligned within the band (forward direction)
//
// mbegin: the state index where the alignment begins (1-indexed)
// seq: sequence to be aligned with the 'suffix' of the model (0-indexed) 
void SeqAlignExtendP::AlignInitFW(
    int mbegin, std::string &seq, AlnResultStruct &result
) {
  assert(mbegin >= 1 && mbegin <= model_.GetProfileLen());
  assert(!seq.empty());
  // calculate null score accociate with the branch
  CalNull(0, seq);
  int i, j;
  // initialize the mmx, imx, and dmx tables
  // note that HMM-style initializetion is different from SW-style
  // also note that the MODEL index is 1-based and the sequence string is 0-based
  mmx_[0][0] = 0; // this is the only thing needs to be done for match-matrix
  // this score is, in fact, never used by the main program. setting the value is
  // only for convinience reason. Note that if mbegin == 1 than B->I0 is used implicitly 
  // instead of M->I
  
  imx_[0][0] = model_.TranScore(mbegin - 1, MI);  
  for(j = 1; j <= Min(band_size_, (int) seq.length()); ++ j) {
    // note that transition score for M->I is added in cell [0][0]
    imx_[0][j] = imx_[0][j - 1]
        + model_.TranScore(mbegin - 1, II) // transition score I->I
        + model_.IEmitScore(mbegin - 1, seq[j - 1]);  // emit score for the amino acids
  }
  // set for the same reason as for imx_[0][0]. Note that if mbegin == 1 then 
  // B->I0 is used implicitly instead of M->D
  dmx_[0][0] = model_.TranScore(mbegin - 1, MD);   
  for(i = 1; i <= Min(band_size_, model_.GetProfileLen() - mbegin + 1); ++ i) {
    dmx_[i][0] = dmx_[i - 1][0]
        + model_.TranScore(mbegin + i - 1, DD);
    // note that there is no emission score for deletion states
  }
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  // fill-up the DP matrix, within the banded regions only
  int m_row = pflen_ - mbegin + 1;      // maximum ending index of rows
  int m_col = seq.length();             // maximum ending index of columns
  int n_row = Min(m_row, m_col + band_size_); // actual ending index of rows, bounded by band size
  int n_col = Min(m_col, m_row + band_size_); // actual ending index of columns, bounded by band size
  SCORETYPE opt_score = INF;
  for(i = 1; i <= n_row; ++ i) {
    pair<int, int> range = GetBandRegionRow(i);
    range.first = range.first < 1 ? 1 : range.first;
    //cout << "range: " << range.first << " " << range.second << endl;
    for(j = range.first; j <= Min(range.second, n_col); ++ j) {
      // note that the score used in HMMER3 is negative log of the probabilities
      // in this case we need to MINIMIZE the sum of the scores
      // Compute matching state
      FillScoreCoreFW(i, j, mbegin, seq, 0);
      // compute the best score here, since no penality is applied for locality
      // insertion and deletion are impossible to be optimal
      if(mmx_[i][j] - null_score_[j - 1] < opt_score)  {
        opt_score = mmx_[i][j] - null_score_[j - 1];
        result.m_score = mmx_[i][j];
        result.n_score = null_score_[j - 1];
        result.m_index = i;
        result.s_index = j;
      }
    }
  }
  result.m_filled = n_row;
  result.s_filled = n_col;
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  return;
}

// EXTEND EXISTING ALIGNMENT WITH NEW SEQUENCE "SEQ"
// mbegin: the state index where the inital alignment begins (1-based)
//         this should be directly adjacent to the identified seed
// m_filled, s_filled: the index describing the region that has already been filled,
//         alignment should be resumed with new sequence by taking the already computed
//         scores in the filled area
// seq: the appended sequence for the new alignment
void SeqAlignExtendP::AlignExtdFW(
    int mbegin,    
    int m_filled, int s_filled,
    std::string &seq,
    AlnResultStruct &result
) {
  assert(mbegin >= 1 && mbegin <= model_.GetProfileLen());
  assert(!seq.empty());
  // calculate null score accociate with the branch
  CalNull(s_filled, seq);
  int i, j;
  // calculate the boundary for the new alignment extension
  int m_row = pflen_ - mbegin + 1;      // maximum ending index of rows
  int m_col = s_filled + seq.length();  // maximum ending index of columns
  int n_row = Min(m_row, m_col + band_size_); // actual ending index of rows, bounded by band size
  int n_col = Min(m_col, m_row + band_size_); // actual ending index of columns, bounded by band size
  assert(n_row >= m_filled && n_col >= s_filled);
  SCORETYPE opt_score = INF;
  // determine the begin of the row need to be filled
  int b = s_filled - band_size_;
  b = b < 1 ? 1 : b;
  for(i = b; i <= n_row; ++ i) {
    // determine the edge begin index
    pair<int, int> range = GetBandRegionRow(i);
    range.first = range.first > s_filled + 1 ? range.first : s_filled + 1;
    range.second = range.second > n_col ? n_col : range.second;
    for(j = range.first; j <= range.second; ++ j) {
      // fill in the matrix
      FillScoreCoreFW(i, j, mbegin, seq, s_filled);
      // compute the best score here, since no penality is applied for locality
      // insertion and deletion are impossible to be optimal
      if(mmx_[i][j] - null_score_[j - 1] < opt_score)  {
        opt_score = mmx_[i][j] - null_score_[j - 1];
        result.m_score = mmx_[i][j];
        result.n_score = null_score_[j - 1];
        result.m_index = i;
        result.s_index = j;
      }
    }
  }
  result.m_filled = n_row;
  result.s_filled = n_col;
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < s_filled + (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  return;
}

// the very first alignment with a sequence intput, the entire
// sequence is to be aligned within the band (reverse direction)
//
// mbegin: the state index where the alignment begins (1-indexed)
// seq: sequence to be aligned with the 'prefix' of the model (0-indexed) 
//      the sequence should be reversed
void SeqAlignExtendP::AlignInitRE(
    int mbegin, std::string &seq, AlnResultStruct &result
) {
  assert(mbegin >= 1 && mbegin <= model_.GetProfileLen());
  assert(!seq.empty());
  // calculate null score accociate with the branch
  CalNull(0, seq);
  int i, j;
  // initialize the mmx, imx, and dmx tables
  // note that HMM-style initializetion is different from SW-style
  // also note that the MODEL index is 1-based and the sequence string is 0-based
  mmx_[0][0] = 0; // this is the only thing needs to be done for match-matrix
  // this score is, in fact, never used by the main program. setting the value is
  // only for convinience reason. Note that if mbegin == 1 than B->I0 is used implicitly 
  // instead of M->I
  //cout << "AlignInitRE: begin initializing scoring matrix" << endl;
  imx_[0][0] = model_.TranScore(mbegin, IM);  
  for(j = 1; j <= Min(band_size_, (int) seq.length()); ++ j) {
    // note that transition score for M->I is added in cell [0][0]
    imx_[0][j] = imx_[0][j - 1]
        + model_.TranScore(mbegin, II) // transition score I->I
        + model_.IEmitScore(mbegin, seq[j - 1]);  // emit score for the amino acids
  }
  // set for the same reason as for imx_[0][0]. Note that if mbegin == 1 then 
  // B->I0 is used implicitly instead of M->D
  dmx_[0][0] = model_.TranScore(mbegin, DM);   
  for(i = 1; i <= Min(band_size_, mbegin); ++ i) {
    dmx_[i][0] = dmx_[i - 1][0]
        + model_.TranScore(mbegin - i, DD);
    // note that there is no emission score for deletion states
  }
  //cout << "AlignInitRE: finished initializing scoring matrix" << endl;
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  // fill-up the DP matrix, within the banded regions only
  int m_row = mbegin;                   // maximum ending index of rows
  int m_col = seq.length();             // maximum ending index of columns
  int n_row = Min(m_row, m_col + band_size_); // actual ending index of rows, bounded by band size
  int n_col = Min(m_col, m_row + band_size_); // actual ending index of columns, bounded by band size
  SCORETYPE opt_score = INF;
  for(i = 1; i <= n_row; ++ i) {
    pair<int, int> range = GetBandRegionRow(i);
    range.first = range.first < 1 ? 1 : range.first;
    //cout << "range: " << range.first << " " << range.second << endl;
    for(j = range.first; j <= Min(range.second, n_col); ++ j) {
      // note that the score used in HMMER3 is negative log of the probabilities
      // in this case we need to MINIMIZE the sum of the scores
      //cout << i << "  " << j << endl;
      // Compute matching state
      FillScoreCoreRE(i, j, mbegin, seq, 0);
      // compute the best score here, since no penality is applied for locality
      // insertion and deletion are impossible to be optimal
      if(mmx_[i][j] - null_score_[j - 1] < opt_score)  {
        opt_score = mmx_[i][j] - null_score_[j - 1];
        result.m_score = mmx_[i][j];
        result.n_score = null_score_[j - 1];
        result.m_index = i;
        result.s_index = j;
      }
    }
  }
  result.m_filled = n_row;
  result.s_filled = n_col;
  //cout << "AlignInitRE dimension: " << n_row << " " << n_col << endl;
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  return;
}

// EXTEND EXISTING ALIGNMENT WITH NEW SEQUENCE "SEQ" (reverse direction)
// mbegin: the state index where the inital alignment begins (1-based)
//         this should be directly adjacent to the identified seed
// m_filled, s_filled: the index describing the region that has already been filled,
//         alignment should be resumed with new sequence by taking the already computed
//         scores in the filled area
// seq: the appended sequence for the new alignment
void SeqAlignExtendP::AlignExtdRE(
    int mbegin,    
    int m_filled, int s_filled,
    std::string &seq,
    AlnResultStruct &result
) {
  assert(mbegin >= 1 && mbegin <= model_.GetProfileLen());
  assert(!seq.empty());
  // calculate null score accociate with the branch
  CalNull(s_filled, seq);
  int i, j;
  // calculate the boundary for the new alignment extension
  int m_row = mbegin + 1;               // maximum ending index of rows
  int m_col = s_filled + seq.length();  // maximum ending index of columns
  int n_row = Min(m_row, m_col + band_size_); // actual ending index of rows, bounded by band size
  int n_col = Min(m_col, m_row + band_size_); // actual ending index of columns, bounded by band size
  //cout << "AlignExtdRE extension dimensions:  " << n_row << " " << n_col << " " << m_filled << "  " << s_filled << endl;
  assert(n_row >= m_filled && n_col >= s_filled);
  SCORETYPE opt_score = INF;
  // determine the begin of the row need to be filled
  int b = s_filled - band_size_;
  b = b < 1 ? 1 : b;
  for(i = b; i <= n_row; ++ i) {
    // determine the edge begin index
    pair<int, int> range = GetBandRegionRow(i);
    range.first = range.first > s_filled + 1 ? range.first : s_filled + 1;
    range.second = range.second > n_col ? n_col : range.second;
    for(j = range.first; j <= range.second; ++ j) {
      // fill in the matrix
      FillScoreCoreRE(i, j, mbegin, seq, s_filled);
      // compute the best score here, since no penality is applied for locality
      // insertion and deletion are impossible to be optimal
      if(mmx_[i][j] - null_score_[j - 1] < opt_score)  {
        opt_score = mmx_[i][j] - null_score_[j - 1];
        result.m_score = mmx_[i][j];
        result.n_score = null_score_[j - 1];
        result.m_index = i;
        result.s_index = j;
      }
    }
  }
  result.m_filled = n_row;
  result.s_filled = n_col;
  // DEBUG
  /*
  cout << "NULL:" << endl;
  for(i = 0; i < s_filled + (int) seq.length(); ++ i) {
    cout << null_score_[i] << " ";
  }
  cout << endl;
  cout << "MMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << mmx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "IMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << imx_[i][j] << " ";
    }
    cout << endl;
  }
  cout << "DMX:" << endl;
  for(i = 0; i < pflen_ + 1; ++ i) {
    for(j = 0; j < pflen_ + band_size_ + 1; ++ j) {
      cout << dmx_[i][j] << " ";
    }
    cout << endl;
  }
  */
  // END DEBUG
  return;
}

/*
SCORETYPE SeqAlignExtendP::Align(
    enum DIRECTION direction, 
    HMMProfile &model, int m_begin,  
    std::string &seq, int s_begin,
    
)  {

}
*/
