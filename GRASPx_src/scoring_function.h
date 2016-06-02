#include <assert.h>
#include <iostream>
#include <map>
#include <cmath>
#include <cstdlib>
#include <vector>

#ifndef _SCORING_FUNCTION_H_
#define _SCORING_FUNCTION_H_

#ifndef _MAX_VALUE
#define _MAX_VALUE 30000
#endif

enum StringType {DNA, RNA, PROTEIN};
enum MatrixName {BLOSUM62, BLOSUM80, BLOSUM90, BLOSUM50, BLOSUM45, PAM250, PAM70, PAM30};

template<typename SCORETYPE>
class ScoringFunction  {
 public:
  ScoringFunction(void);
  ScoringFunction(SCORETYPE in_match, SCORETYPE in_mismatch, SCORETYPE in_gap_extend, SCORETYPE in_gap_open);
  ScoringFunction(enum StringType mol, enum MatrixName matrix, SCORETYPE in_gap_extend, SCORETYPE in_gap_open);
  ~ScoringFunction(void);
  SCORETYPE CheckMatchScore(char p, char q); 
  SCORETYPE GetGapOpen(void);
  SCORETYPE GetGapExtend(void);
  double GetKarlinStatK(void);
  double GetKarlinStatLambda(void);
  double ComputeBitScore(SCORETYPE s);
  double ComputeEValue(int m, int n, SCORETYPE s);
  SCORETYPE ComputeRawScore(long int m, long int n, double e);
  ScoringFunction& operator = (const ScoringFunction &in_obj);
  SCORETYPE GetAveMatch(void) {
    return ave_match;
  }
  SCORETYPE CalMatchScore(std::string& s1, std::string& s2);
 protected:
  bool has_matrix_initialized;
  SCORETYPE match; 
  SCORETYPE mismatch; 
  SCORETYPE gap_extend; 
  SCORETYPE gap_open;
  SCORETYPE ave_match;
  std::vector< std::vector<SCORETYPE> > score_matrix;
  double Karlin_stat_K;
  double Karlin_stat_Lambda;
  std::map<char, int> alphabet_hash;
  void compute_ave_match(void);
  void fill_matrix_BLOSUM62(void);
  void fill_matrix_BLOSUM80(void);
  void fill_matrix_BLOSUM90(void);
  void fill_matrix_BLOSUM50(void);
  void fill_matrix_BLOSUM45(void);
  void fill_matrix_PAM250(void);
  void fill_matrix_PAM70(void);
  void fill_matrix_PAM30(void);
  void get_Karlin_stats_BLOSUM62(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_BLOSUM80(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_BLOSUM90(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_BLOSUM50(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_BLOSUM45(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_PAM250(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_PAM70(SCORETYPE gap_extend, SCORETYPE gap_open);
  void get_Karlin_stats_PAM30(SCORETYPE gap_extend, SCORETYPE gap_open);
  void assign_K_and_Lambda(const std::vector< std::vector<double> > Karlin_matrix, 
    const SCORETYPE gap_extend, const SCORETYPE gap_open);
  
};

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::compute_ave_match(void)  {
  if(!has_matrix_initialized)  {
    ave_match = match;
  } else  {
    SCORETYPE total_score = 0;
    unsigned int i;
    for(i = 0; i < score_matrix.size(); ++ i) {
      total_score += score_matrix[i][i];
    }
    ave_match = (SCORETYPE) (total_score / i);
  }
  return;
}

template<typename SCORETYPE>
SCORETYPE ScoringFunction<SCORETYPE>::CalMatchScore(std::string& s1, std::string& s2)  {
  if(s1.length() != s2.length())  {
    std::cout << "Warnning: ScoringFunction::CalMatchScore incompatible lengths \
      of sequences." << std::endl;
    return 0;
  }
  SCORETYPE score = 0;
  for(unsigned int i = 0; i < s1.length(); ++ i) {
    score += (SCORETYPE) CheckMatchScore(s1[i], s2[i]); 
  }
  return score;
}

template<typename SCORETYPE> SCORETYPE ScoringFunction<SCORETYPE>::CheckMatchScore(char p, char q)  {
  if(!has_matrix_initialized)  {
    return (p == q ? match : mismatch);
  }  else  {
    int i_index, j_index;
    if(alphabet_hash.find(p) == alphabet_hash.end())  {
      // treating the character as 'X'
      i_index = 23;
    } else  {
      i_index = alphabet_hash[p];
    }
    if(alphabet_hash.find(q) == alphabet_hash.end())  {
      // treating the character as 'X'
      j_index = 23;
    } else  {
      j_index = alphabet_hash[q];
    }
    return score_matrix[i_index][j_index];
  }
}

template<typename SCORETYPE>
double ScoringFunction<SCORETYPE>::GetKarlinStatK(void)  {
  assert(has_matrix_initialized);
  return Karlin_stat_K;
}

template<typename SCORETYPE>
double ScoringFunction<SCORETYPE>::GetKarlinStatLambda(void)  {
  assert(has_matrix_initialized);
  return Karlin_stat_Lambda;
}

template<typename SCORETYPE> ScoringFunction<SCORETYPE>::ScoringFunction(void)  {
  // a naive one, for backward compatibility
  
  has_matrix_initialized = false;
  
  match = 3;
  mismatch = -4;
  gap_extend = -2;
  gap_open = -10;
  compute_ave_match();
  return;
}

template<typename SCORETYPE>
ScoringFunction<SCORETYPE>::ScoringFunction(SCORETYPE in_match, SCORETYPE in_mismatch, SCORETYPE in_gap_extend, SCORETYPE in_gap_open)  {
  // initializing simple scoring function

  has_matrix_initialized = false;
  
  match = in_match;
  mismatch = in_mismatch;
  gap_extend = in_gap_extend;
  gap_open = in_gap_open;
  compute_ave_match();
  return;
}

template<typename SCORETYPE>
ScoringFunction<SCORETYPE>::ScoringFunction(enum StringType mol, enum MatrixName matrix, SCORETYPE in_gap_extend, SCORETYPE in_gap_open)  {
  
  has_matrix_initialized = true;
  
  gap_extend = in_gap_extend;
  gap_open = in_gap_open;
  
  if(mol == PROTEIN)  {
    // the amino acids order is defined as the scoring matrix used in BLAST
    alphabet_hash = {
      {'A', 0},  {'R', 1},  {'N', 2},  {'D', 3},  {'C', 4},  {'Q', 5},  {'E', 6},
      {'G', 7},  {'H', 8},  {'I', 9},  {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13},
      {'P', 14}, {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19}, {'B', 20},
      {'J', 21}, {'Z', 22}, {'X', 23}, {'*', 24}
    };
    score_matrix.resize(25);
    for(unsigned int i = 0; i < 25; ++ i)  {
      score_matrix[i].resize(25);
    }
    switch(matrix)  {
      case BLOSUM62:
        fill_matrix_BLOSUM62();
        get_Karlin_stats_BLOSUM62(in_gap_extend, in_gap_open);
        break;
      case BLOSUM80:
        fill_matrix_BLOSUM80();
        get_Karlin_stats_BLOSUM80(in_gap_extend, in_gap_open);
        break;
      case BLOSUM90:
        fill_matrix_BLOSUM90();
        get_Karlin_stats_BLOSUM90(in_gap_extend, in_gap_open);
        break;
      case BLOSUM50:
        fill_matrix_BLOSUM50();
        get_Karlin_stats_BLOSUM50(in_gap_extend, in_gap_open);
        break;
      case BLOSUM45:
        fill_matrix_BLOSUM45();
        get_Karlin_stats_BLOSUM45(in_gap_extend, in_gap_open);
        break;
      case PAM250:
        fill_matrix_PAM250();
        get_Karlin_stats_PAM250(in_gap_extend, in_gap_open);
        break;
      case PAM70:
        fill_matrix_PAM70();
        get_Karlin_stats_PAM70(in_gap_extend, in_gap_open);
        break;
      case PAM30:
        fill_matrix_PAM30();
        get_Karlin_stats_PAM30(in_gap_extend, in_gap_open);
        break;
      default:
        std::cout << "No such scoring function is currently supported." << std::endl;
        break;
    }
  }  else if(mol == DNA)  {
    alphabet_hash = {
      {'A', 0},  {'C', 0},  {'G', 0},  {'T', 3},  {'N', 4}
    };
    std::cout << "Currently DNA alignment scoring matrix is not supported." << std::endl;
    exit(0);
  }  else if(mol == RNA)  {
    alphabet_hash = {
      {'A', 0},  {'C', 0},  {'G', 0},  {'U', 3},  {'N', 4}
    };
    std::cout << "Currently RNA alignment scoring matrix is not supported." << std::endl;
    exit(0);
  }  else  {
    std::cout << "Currently generic string alignment scoring matrix is not supported." << std::endl;
    exit(0);
  }
  compute_ave_match();
  return;
}

template<typename SCORETYPE> ScoringFunction<SCORETYPE>::~ScoringFunction(void)  {
  return;
}

template<typename SCORETYPE>
double ScoringFunction<SCORETYPE>::ComputeBitScore(SCORETYPE s)  {
  assert(has_matrix_initialized);
  //std::cout << Karlin_stat_Lambda << " " << Karlin_stat_K << std::endl;
  double bit_s = (Karlin_stat_Lambda * (double) s - log(Karlin_stat_K)) / log(2);
  return bit_s;
}

template<typename SCORETYPE>
double ScoringFunction<SCORETYPE>::ComputeEValue(int m, int n, SCORETYPE s) {
  assert(has_matrix_initialized);
  assert(m > 0 && n > 0);
  double e_value = Karlin_stat_K * m * n * exp(-Karlin_stat_Lambda * s);
  return e_value;
}

template<typename SCORETYPE>
SCORETYPE ScoringFunction<SCORETYPE>::ComputeRawScore(long int m, long int n, double e) {
  assert(has_matrix_initialized);
  assert(m > 0 && n > 0);
  SCORETYPE raw_s = static_cast<SCORETYPE>((log(Karlin_stat_K) + log(m) + log (n) - log(e)) / Karlin_stat_Lambda);
  return raw_s;
}

template<typename SCORETYPE>
SCORETYPE ScoringFunction<SCORETYPE>::GetGapOpen()  {
  return gap_open;
}

template<typename SCORETYPE>
SCORETYPE ScoringFunction<SCORETYPE>::GetGapExtend()  {
  return gap_extend;
}

template<typename SCORETYPE> ScoringFunction<SCORETYPE>& ScoringFunction<SCORETYPE>::operator = (const ScoringFunction &in_obj)  {
  this->match = in_obj.match;
  this->mismatch = in_obj.match;
  this->gap_extend = in_obj.gap_extend;
  this->gap_open = in_obj.gap_open;
  return *this;
}


template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_BLOSUM62(void)  {
  score_matrix = {
    { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4},
    {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4},
    {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4},
    {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4},
    { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4},
    {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4},
    {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4},
    { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4},
    {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4},
    {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4},
    {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4},
    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4},
    {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4},
    { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4},
    {-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4},
    {-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4},
    {-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4},
    {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_BLOSUM80(void)  {
  score_matrix = {
    { 5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -6},
    {-2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3,  0, -1, -6},
    {-2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4,  5, -4,  0, -1, -6},
    {-2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4,  5, -5,  1, -1, -6},
    {-1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1, -6},
    {-1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3,  0, -3,  4, -1, -6},
    {-1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3,  1, -4,  5, -1, -6},
    { 0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6},
    {-2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -1, -4,  0, -1, -6},
    {-2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -4,  3, -4, -1, -6},
    {-2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -4,  3, -3, -1, -6},
    {-1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3,  1, -1, -6},
    {-1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -3,  2, -1, -1, -6},
    {-3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -4,  0, -4, -1, -6},
    {-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -2, -4, -2, -1, -6},
    { 1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0, -3,  0, -1, -6},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -1, -1, -1, -1, -6},
    {-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5, -3, -3, -1, -6},
    {-2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -3, -2, -3, -1, -6},
    { 0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -4,  2, -3, -1, -6},
    {-2, -1,  5,  5, -4,  0,  1, -1, -1, -4, -4, -1, -3, -4, -2,  0, -1, -5, -3, -4,  5, -4,  0, -1, -6},
    {-2, -3, -4, -5, -2, -3, -4, -5, -4,  3,  3, -3,  2,  0, -4, -3, -1, -3, -2,  2, -4,  3, -3, -1, -6},
    {-1,  0,  0,  1, -4,  4,  5, -3,  0, -4, -3,  1, -1, -4, -2,  0, -1, -3, -3, -3,  0, -3,  5, -1, -6},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6},
    {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_BLOSUM90(void)  {
  score_matrix = {
    { 5, -2, -2, -3, -1, -1, -1,  0, -2, -2, -2, -1, -2, -3, -1,  1,  0, -4, -3, -1, -2, -2, -1, -1, -6},
    {-2,  6, -1, -3, -5,  1, -1, -3,  0, -4, -3,  2, -2, -4, -3, -1, -2, -4, -3, -3, -2, -3,  0, -1, -6},
    {-2, -1,  7,  1, -4,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -5, -3, -4,  5, -4, -1, -1, -6},
    {-3, -3,  1,  7, -5, -1,  1, -2, -2, -5, -5, -1, -4, -5, -3, -1, -2, -6, -4, -5,  5, -5,  1, -1, -6},
    {-1, -5, -4, -5,  9, -4, -6, -4, -5, -2, -2, -4, -2, -3, -4, -2, -2, -4, -4, -2, -4, -2, -5, -1, -6},
    {-1,  1,  0, -1, -4,  7,  2, -3,  1, -4, -3,  1,  0, -4, -2, -1, -1, -3, -3, -3, -1, -3,  5, -1, -6},
    {-1, -1, -1,  1, -6,  2,  6, -3, -1, -4, -4,  0, -3, -5, -2, -1, -1, -5, -4, -3,  1, -4,  5, -1, -6},
    { 0, -3, -1, -2, -4, -3, -3,  6, -3, -5, -5, -2, -4, -5, -3, -1, -3, -4, -5, -5, -2, -5, -3, -1, -6},
    {-2,  0,  0, -2, -5,  1, -1, -3,  8, -4, -4, -1, -3, -2, -3, -2, -2, -3,  1, -4, -1, -4,  0, -1, -6},
    {-2, -4, -4, -5, -2, -4, -4, -5, -4,  5,  1, -4,  1, -1, -4, -3, -1, -4, -2,  3, -5,  3, -4, -1, -6},
    {-2, -3, -4, -5, -2, -3, -4, -5, -4,  1,  5, -3,  2,  0, -4, -3, -2, -3, -2,  0, -5,  4, -4, -1, -6},
    {-1,  2,  0, -1, -4,  1,  0, -2, -1, -4, -3,  6, -2, -4, -2, -1, -1, -5, -3, -3, -1, -3,  1, -1, -6},
    {-2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7, -1, -3, -2, -1, -2, -2,  0, -4,  2, -2, -1, -6},
    {-3, -4, -4, -5, -3, -4, -5, -5, -2, -1,  0, -4, -1,  7, -4, -3, -3,  0,  3, -2, -4,  0, -4, -1, -6},
    {-1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8, -2, -2, -5, -4, -3, -3, -4, -2, -1, -6},
    { 1, -1,  0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2,  5,  1, -4, -3, -2,  0, -3, -1, -1, -6},
    { 0, -2,  0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2,  1,  6, -4, -2, -1, -1, -2, -1, -1, -6},
    {-4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2,  0, -5, -4, -4, 11,  2, -3, -6, -3, -4, -1, -6},
    {-3, -3, -3, -4, -4, -3, -4, -5,  1, -2, -2, -3, -2,  3, -4, -3, -2,  2,  8, -3, -4, -2, -3, -1, -6},
    {-1, -3, -4, -5, -2, -3, -3, -5, -4,  3,  0, -3,  0, -2, -3, -2, -1, -3, -3,  5, -4,  1, -3, -1, -6},
    {-2, -2,  5,  5, -4, -1,  1, -2, -1, -5, -5, -1, -4, -4, -3,  0, -1, -6, -4, -4,  5, -5,  0, -1, -6},
    {-2, -3, -4, -5, -2, -3, -4, -5, -4,  3,  4, -3,  2,  0, -4, -3, -2, -3, -2,  1, -5,  4, -4, -1, -6},
    {-1,  0, -1,  1, -5,  5,  5, -3,  0, -4, -4,  1, -2, -4, -2, -1, -1, -4, -3, -3,  0, -4,  5, -1, -6},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6},
    {-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_BLOSUM50(void)  {
  score_matrix = {
    { 5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -5},
    {-2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1, -3,  0, -1, -5},
    {-1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5, -4,  0, -1, -5},
    {-2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6, -4,  1, -1, -5},
    {-1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -2, -3, -1, -5},
    {-1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0, -3,  4, -1, -5},
    {-1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1, -3,  5, -1, -5},
    { 0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -4, -2, -1, -5},
    {-2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0, -3,  0, -1, -5},
    {-1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4,  4, -3, -1, -5},
    {-2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4,  4, -3, -1, -5},
    {-1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0, -3,  1, -1, -5},
    {-1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3,  2, -1, -1, -5},
    {-3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4,  1, -4, -1, -5},
    {-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -3, -1, -1, -5},
    { 1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0, -3,  0, -1, -5},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -1, -5},
    {-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -2, -1, -5},
    {-2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -1, -2, -1, -5},
    { 0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3,  2, -3, -1, -5},
    {-2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6, -4,  1, -1, -5},
    {-2, -3, -4, -4, -2, -3, -3, -4, -3,  4,  4, -3,  2,  1, -3, -3, -1, -2, -1,  2, -4,  4, -3, -1, -5},
    {-1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1, -3,  5, -1, -5},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
    {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_BLOSUM45(void)  {
  score_matrix = {
    { 5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1, -1, -1, -5},
    {-2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1, -3,  1, -1, -5},
    {-1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  5, -3,  0, -1, -5},
    {-2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  6, -3,  1, -1, -5},
    {-1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -2, -3, -1, -5},
    {-1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0, -2,  4, -1, -5},
    {-1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1, -3,  5, -1, -5},
    { 0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -5},
    {-2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0, -2,  0, -1, -5},
    {-1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3,  4, -3, -1, -5},
    {-1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3,  4, -2, -1, -5},
    {-1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0, -3,  1, -1, -5},
    {-1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2,  2, -1, -1, -5},
    {-2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3,  1, -3, -1, -5},
    {-1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -3, -1, -1, -5},
    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0, -2,  0, -1, -5},
    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1, -1, -1, -5},
    {-2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -1, -5},
    {-2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2,  0, -2, -1, -5},
    { 0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3,  2, -3, -1, -5},
    {-1, -1,  5,  6, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  5, -3,  1, -1, -5},
    {-1, -3, -3, -3, -2, -2, -3, -4, -2,  4,  4, -3,  2,  1, -3, -2, -1, -2,  0,  2, -3,  4, -2, -1, -5},
    {-1,  1,  0,  1, -3,  4,  5, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  1, -2,  5, -1, -5},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5},
    {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_PAM250(void)  {
  score_matrix = {
    { 2, -2,  0,  0, -2,  0,  0,  1, -1, -1, -2, -1, -1, -3,  1,  1,  1, -6, -3,  0,  0, -1,  0, -1, -8},
    {-2,  6,  0, -1, -4,  1, -1, -3,  2, -2, -3,  3,  0, -4,  0,  0, -1,  2, -4, -2, -1, -3,  0, -1, -8},
    { 0,  0,  2,  2, -4,  1,  1,  0,  2, -2, -3,  1, -2, -3,  0,  1,  0, -4, -2, -2,  2, -3,  1, -1, -8},
    { 0, -1,  2,  4, -5,  2,  3,  1,  1, -2, -4,  0, -3, -6, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8},
    {-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3,  0, -2, -8,  0, -2, -4, -5, -5, -1, -8},
    { 0,  1,  1,  2, -5,  4,  2, -1,  3, -2, -2,  1, -1, -5,  0, -1, -1, -5, -4, -2,  1, -2,  3, -1, -8},
    { 0, -1,  1,  3, -5,  2,  4,  0,  1, -2, -3,  0, -2, -5, -1,  0,  0, -7, -4, -2,  3, -3,  3, -1, -8},
    { 1, -3,  0,  1, -3, -1,  0,  5, -2, -3, -4, -2, -3, -5,  0,  1,  0, -7, -5, -1,  0, -4,  0, -1, -8},
    {-1,  2,  2,  1, -3,  3,  1, -2,  6, -2, -2,  0, -2, -2,  0, -1, -1, -3,  0, -2,  1, -2,  2, -1, -8},
    {-1, -2, -2, -2, -2, -2, -2, -3, -2,  5,  2, -2,  2,  1, -2, -1,  0, -5, -1,  4, -2,  3, -2, -1, -8},
    {-2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6, -3,  4,  2, -3, -3, -2, -2, -1,  2, -3,  5, -3, -1, -8},
    {-1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5,  0, -5, -1,  0,  0, -3, -4, -2,  1, -3,  0, -1, -8},
    {-1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6,  0, -2, -2, -1, -4, -2,  2, -2,  3, -2, -1, -8},
    {-3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9, -5, -3, -3,  0,  7, -1, -4,  2, -5, -1, -8},
    { 1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6,  1,  0, -6, -5, -1, -1, -2,  0, -1, -8},
    { 1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2,  1, -2, -3, -1,  0, -2,  0, -1, -8},
    { 1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3, -5, -3,  0,  0, -1, -1, -1, -8},
    {-6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17,  0, -6, -5, -3, -6, -1, -8},
    {-3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10, -2, -3, -1, -4, -1, -8},
    { 0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4, -2,  2, -2, -1, -8},
    { 0, -1,  2,  3, -4,  1,  3,  0,  1, -2, -3,  1, -2, -4, -1,  0,  0, -5, -3, -2,  3, -3,  2, -1, -8},
    {-1, -3, -3, -3, -5, -2, -3, -4, -2,  3,  5, -3,  3,  2, -2, -2, -1, -3, -1,  2, -3,  5, -2, -1, -8},
    { 0,  0,  1,  3, -5,  3,  3,  0,  2, -2, -3,  0, -2, -5,  0,  0, -1, -6, -4, -2,  2, -2,  3, -1, -8},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -8},
    {-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_PAM70(void)  {
  score_matrix = {
    {  5, -4, -2, -1, -4, -2, -1,  0, -4, -2, -4, -4, -3, -6,  0,  1,  1, -9, -5, -1, -1, -3, -1, -1,-11},
    { -4,  8, -3, -6, -5,  0, -5, -6,  0, -3, -6,  2, -2, -7, -2, -1, -4,  0, -7, -5, -4, -5, -2, -1,-11},
    { -2, -3,  6,  3, -7, -1,  0, -1,  1, -3, -5,  0, -5, -6, -3,  1,  0, -6, -3, -5,  5, -4, -1, -1,-11},
    { -1, -6,  3,  6, -9,  0,  3, -1, -1, -5, -8, -2, -7,-10, -4, -1, -2,-10, -7, -5,  5, -7,  2, -1,-11},
    { -4, -5, -7, -9,  9, -9, -9, -6, -5, -4,-10, -9, -9, -8, -5, -1, -5,-11, -2, -4, -8, -7, -9, -1,-11},
    { -2,  0, -1,  0, -9,  7,  2, -4,  2, -5, -3, -1, -2, -9, -1, -3, -3, -8, -8, -4, -1, -3,  5, -1,-11},
    { -1, -5,  0,  3, -9,  2,  6, -2, -2, -4, -6, -2, -4, -9, -3, -2, -3,-11, -6, -4,  2, -5,  5, -1,-11},
    {  0, -6, -1, -1, -6, -4, -2,  6, -6, -6, -7, -5, -6, -7, -3,  0, -3,-10, -9, -3, -1, -7, -3, -1,-11},
    { -4,  0,  1, -1, -5,  2, -2, -6,  8, -6, -4, -3, -6, -4, -2, -3, -4, -5, -1, -4,  0, -4,  1, -1,-11},
    { -2, -3, -3, -5, -4, -5, -4, -6, -6,  7,  1, -4,  1,  0, -5, -4, -1, -9, -4,  3, -4,  4, -4, -1,-11},
    { -4, -6, -5, -8,-10, -3, -6, -7, -4,  1,  6, -5,  2, -1, -5, -6, -4, -4, -4,  0, -6,  5, -4, -1,-11},
    { -4,  2,  0, -2, -9, -1, -2, -5, -3, -4, -5,  6,  0, -9, -4, -2, -1, -7, -7, -6, -1, -5, -2, -1,-11},
    { -3, -2, -5, -7, -9, -2, -4, -6, -6,  1,  2,  0, 10, -2, -5, -3, -2, -8, -7,  0, -6,  2, -3, -1,-11},
    { -6, -7, -6,-10, -8, -9, -9, -7, -4,  0, -1, -9, -2,  8, -7, -4, -6, -2,  4, -5, -7, -1, -9, -1,-11},
    {  0, -2, -3, -4, -5, -1, -3, -3, -2, -5, -5, -4, -5, -7,  7,  0, -2, -9, -9, -3, -4, -5, -2, -1,-11},
    {  1, -1,  1, -1, -1, -3, -2,  0, -3, -4, -6, -2, -3, -4,  0,  5,  2, -3, -5, -3,  0, -5, -2, -1,-11},
    {  1, -4,  0, -2, -5, -3, -3, -3, -4, -1, -4, -1, -2, -6, -2,  2,  6, -8, -4, -1, -1, -3, -3, -1,-11},
    { -9,  0, -6,-10,-11, -8,-11,-10, -5, -9, -4, -7, -8, -2, -9, -3, -8, 13, -3,-10, -7, -5,-10, -1,-11},
    { -5, -7, -3, -7, -2, -8, -6, -9, -1, -4, -4, -7, -7,  4, -9, -5, -4, -3,  9, -5, -4, -4, -7, -1,-11},
    { -1, -5, -5, -5, -4, -4, -4, -3, -4,  3,  0, -6,  0, -5, -3, -3, -1,-10, -5,  6, -5,  1, -4, -1,-11},
    { -1, -4,  5,  5, -8, -1,  2, -1,  0, -4, -6, -1, -6, -7, -4,  0, -1, -7, -4, -5,  5, -5,  1, -1,-11},
    { -3, -5, -4, -7, -7, -3, -5, -7, -4,  4,  5, -5,  2, -1, -5, -5, -3, -5, -4,  1, -5,  5, -4, -1,-11},
    { -1, -2, -1,  2, -9,  5,  5, -3,  1, -4, -4, -2, -3, -9, -2, -2, -3,-10, -7, -4,  1, -4,  5, -1,-11},
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-11},
    {-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,-11,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::fill_matrix_PAM30(void)  {
  score_matrix = {
    {  6, -7, -4, -3, -6, -4, -2, -2, -7, -5, -6, -7, -5, -8, -2,  0, -1,-13, -8, -2, -3, -6, -3, -1,-17},
    { -7,  8, -6,-10, -8, -2, -9, -9, -2, -5, -8,  0, -4, -9, -4, -3, -6, -2,-10, -8, -7, -7, -4, -1,-17},
    { -4, -6,  8,  2,-11, -3, -2, -3,  0, -5, -7, -1, -9, -9, -6,  0, -2, -8, -4, -8,  6, -6, -3, -1,-17},
    { -3,-10,  2,  8,-14, -2,  2, -3, -4, -7,-12, -4,-11,-15, -8, -4, -5,-15,-11, -8,  6,-10,  1, -1,-17},
    { -6, -8,-11,-14, 10,-14,-14, -9, -7, -6,-15,-14,-13,-13, -8, -3, -8,-15, -4, -6,-12, -9,-14, -1,-17},
    { -4, -2, -3, -2,-14,  8,  1, -7,  1, -8, -5, -3, -4,-13, -3, -5, -5,-13,-12, -7, -3, -5,  6, -1,-17},
    { -2, -9, -2,  2,-14,  1,  8, -4, -5, -5, -9, -4, -7,-14, -5, -4, -6,-17, -8, -6,  1, -7,  6, -1,-17},
    { -2, -9, -3, -3, -9, -7, -4,  6, -9,-11,-10, -7, -8, -9, -6, -2, -6,-15,-14, -5, -3,-10, -5, -1,-17},
    { -7, -2,  0, -4, -7,  1, -5, -9,  9, -9, -6, -6,-10, -6, -4, -6, -7, -7, -3, -6, -1, -7, -1, -1,-17},
    { -5, -5, -5, -7, -6, -8, -5,-11, -9,  8, -1, -6, -1, -2, -8, -7, -2,-14, -6,  2, -6,  5, -6, -1,-17},
    { -6, -8, -7,-12,-15, -5, -9,-10, -6, -1,  7, -8,  1, -3, -7, -8, -7, -6, -7, -2, -9,  6, -7, -1,-17},
    { -7,  0, -1, -4,-14, -3, -4, -7, -6, -6, -8,  7, -2,-14, -6, -4, -3,-12, -9, -9, -2, -7, -4, -1,-17},
    { -5, -4, -9,-11,-13, -4, -7, -8,-10, -1,  1, -2, 11, -4, -8, -5, -4,-13,-11, -1,-10,  0, -5, -1,-17},
    { -8, -9, -9,-15,-13,-13,-14, -9, -6, -2, -3,-14, -4,  9,-10, -6, -9, -4,  2, -8,-10, -2,-13, -1,-17},
    { -2, -4, -6, -8, -8, -3, -5, -6, -4, -8, -7, -6, -8,-10,  8, -2, -4,-14,-13, -6, -7, -7, -4, -1,-17},
    {  0, -3,  0, -4, -3, -5, -4, -2, -6, -7, -8, -4, -5, -6, -2,  6,  0, -5, -7, -6, -1, -8, -5, -1,-17},
    { -1, -6, -2, -5, -8, -5, -6, -6, -7, -2, -7, -3, -4, -9, -4,  0,  7,-13, -6, -3, -3, -5, -6, -1,-17},
    {-13, -2, -8,-15,-15,-13,-17,-15, -7,-14, -6,-12,-13, -4,-14, -5,-13, 13, -5,-15,-10, -7,-14, -1,-17},
    { -8,-10, -4,-11, -4,-12, -8,-14, -3, -6, -7, -9,-11,  2,-13, -7, -6, -5, 10, -7, -6, -7, -9, -1,-17},
    { -2, -8, -8, -8, -6, -7, -6, -5, -6,  2, -2, -9, -1, -8, -6, -6, -3,-15, -7,  7, -8,  0, -6, -1,-17},
    { -3, -7,  6,  6,-12, -3,  1, -3, -1, -6, -9, -2,-10,-10, -7, -1, -3,-10, -6, -8,  6, -8,  0, -1,-17},
    { -6, -7, -6,-10, -9, -5, -7,-10, -7,  5,  6, -7,  0, -2, -7, -8, -5, -7, -7,  0, -8,  6, -6, -1,-17},
    { -3, -4, -3,  1,-14,  6,  6, -5, -1, -6, -7, -4, -5,-13, -4, -5, -6,-14, -9, -6,  0, -6,  6, -1,-17},
    { -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,-17},
    {-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,  1}
  };
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_BLOSUM62(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 11;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.318,  0.134,  0.401},
    {11,  2,  0.297,  0.082,  0.27},
    {10,  2,  0.291,  0.075,  0.23},
    {12,  1,  0.283,  0.059,  0.19},
    {9,  2,  0.279,  0.058,  0.19},
    {8,  2,  0.264,  0.045,  0.15},
    {11,  1,  0.267,  0.041,  0.14},
    {10,  1,  0.243,  0.024,  0.10},
    {7,  2,  0.239,  0.027,  0.10},
    {6,  2,  0.201,  0.012,  0.061},
    {9,  1,  0.206,  0.010,  0.052}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_BLOSUM80(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 10;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.343,  0.177,  0.657},
    {25,  2,  0.342,  0.17,  0.66},
    {13,  2,  0.336,  0.15,  0.57},
    {9,  2,  0.319,  0.11,  0.42},
    {11,  1,  0.314,  0.095,  0.35},
    {8,  2,  0.308,  0.090,  0.35},
    {10,  1,  0.299,  0.071,  0.27},
    {7,  2,  0.293,  0.070,  0.27},
    {9,  1,  0.279,  0.048,  0.20},
    {6,  2,  0.268,  0.045,  0.19}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_BLOSUM90(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 8;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.335,  0.190,  0.755},
     {9,  2,  0.310,  0.12,  0.46},
     {11,  1,  0.302,  0.093,  0.39},
     {8,  2,  0.300,  0.099,  0.39},
     {7,  2,  0.283,  0.072,  0.30},
     {10,  1,  0.290,  0.075,  0.28},
     {6,  2,  0.259,  0.048,  0.22},
     {9,  1,  0.265,  0.044,  0.20}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_BLOSUM50(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 16;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.232,  0.112,  0.336},
    {16,  2,  0.215,  0.066,  0.20},
    {13,  3,  0.212,  0.063,  0.19},
    {19,  1,  0.212,  0.57,  0.18},
    {15,  2,  0.210,  0.058,  0.17},
    {12,  3,  0.206,  0.055,  0.17},
    {18,  1,  0.207,  0.050,  0.15},
    {14,  2,  0.202,  0.045,  0.14},
    {11,  3,  0.197,  0.042,  0.14},
    {17,  1,  0.198,  0.037,  0.12},
    {13,  2,  0.193,  0.035,  0.12},
    {10,  3,  0.186,  0.031,  0.11},
    {16,  1,  0.186,  0.025,  0.10},
    {12,  2,  0.181,  0.025,  0.095},
    {9,  3,  0.172,  0.022,  0.082},
    {15,  1,  0.171,  0.015,  0.063}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_BLOSUM45(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 14;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.229,  0.092,  0.251},
     {13,  3,  0.207,  0.049,  0.14},
     {16,  2,  0.210,  0.051,  0.14},
     {15,  2,  0.203,  0.041,  0.12},
     {19,  1,  0.205,  0.040,  0.11},
     {12,  3,  0.199,  0.039,  0.11},
     {18,  1,  0.198,  0.032,  0.10},
     {14,  2,  0.195,  0.032,  0.10},
     {11,  3,  0.190,  0.031,  0.095},
     {13,  2,  0.185,  0.024,  0.084},
     {17,  1,  0.189,  0.024,  0.078},
     {10,  3,  0.179,  0.023,  0.075},
     {16,  1,  0.176,  0.016,  0.063},
     {12,  2,  0.171,  0.016,  0.061} 
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_PAM250(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 16;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.218,  0.0877,  0.287},
     {15,  3,  0.205,  0.049,  0.13},
     {17,  2,  0.204,  0.047,  0.12},
     {14,  3,  0.200,  0.043,  0.12},
     {21,  1,  0.204,  0.045,  0.11},
     {16,  2,  0.198,  0.038,  0.11},
     {20,  1,  0.199,  0.037,  0.10},
     {13,  3,  0.194,  0.036,  0.10},
     {15,  2,  0.191,  0.031,  0.087},
     {12,  3,  0.186,  0.029,  0.085},
     {19,  1,  0.192,  0.029,  0.083},
     {14,  2,  0.182,  0.024,  0.073},
     {18,  1,  0.183,  0.021,  0.070},
     {11,  3,  0.174,  0.020,  0.070},
     {13,  2,  0.171,  0.017,  0.059},
     {17,  1,  0.171,  0.014,  0.052}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_PAM30(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 7;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.336,  0.277,  1.82},
     {10,  1,  0.309,  0.15,  0.88},
     {7,  2,  0.305,  0.15,  0.87},
     {6,  2,  0.287,  0.11,  0.68},
     {9,  1,  0.294,  0.11,  0.61},
     {5,  2,  0.264,  0.079,  0.45},
     {8,  1,  0.270,  0.072,  0.40}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::get_Karlin_stats_PAM70(SCORETYPE gap_extend, SCORETYPE gap_open)  {
  
  unsigned int i;
  unsigned int num_accepted_groups = 7;
  
  std::vector< std::vector<double> > precomputed_parameters;
  precomputed_parameters.resize(num_accepted_groups);
  for(i = 0; i < num_accepted_groups; ++ i)  {
    precomputed_parameters.resize(5);
  }
  
  precomputed_parameters = {
    {_MAX_VALUE,  _MAX_VALUE,  0.328,  0.222, 1.11},
     {8,  2,  0.301,  0.12,  0.54},
     {11,  1,  0.305,  0.12,  0.52},
     {7,  2,  0.286,  0.093,  0.43},
     {10,  1,  0.291,  0.91,  0.41},
     {6,  2,  0.264,  0.064,  0.29},
     {9,  1,  0.270,  0.060,  0.28}
  };
  
  assign_K_and_Lambda(precomputed_parameters, gap_extend, gap_open);
  
  return;
}

template<typename SCORETYPE>
void ScoringFunction<SCORETYPE>::assign_K_and_Lambda(
  const std::vector< std::vector<double> > Karlin_matrix, const SCORETYPE gap_extend, const SCORETYPE gap_open
)  {

  unsigned int i;
  for(i = 1; i < Karlin_matrix.size(); ++ i)  {
    if(gap_open == -(SCORETYPE) Karlin_matrix[i][0] && gap_extend == -(SCORETYPE) Karlin_matrix[i][1])  {
      Karlin_stat_K = Karlin_matrix[i][3];
      Karlin_stat_Lambda = Karlin_matrix[i][2];
      return;
    }
  }
  std::cout << "The input gap open and extension penalty is not supported for the selected scoring matrix" << std::endl;
  std::cout << "Allowed gap open and extension penalty combinations are:" << std::endl;
  std::cout << "(gap_open,gap_extension)" << std::endl;
  for(i = 1; i < Karlin_matrix.size(); ++ i)  { 
    std::cout << "(" << -(SCORETYPE) Karlin_matrix[i][0] << "," << -(SCORETYPE) Karlin_matrix[i][1] << ")" << std::endl;
  }
  std::cout << "Please rerun Search by choosing one of the allowed combinations using --gap_open and --gap_extension options." << std::endl;
  exit(0);
}

#endif
