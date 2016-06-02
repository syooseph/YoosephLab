#include <boost/algorithm/string.hpp>
#include "reachable_reads.h"
#include "database_index.h"

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <cstdlib>
#include <map>
#include <cmath>

#ifndef _HMM_PROFILE_
#define _HMM_PROFILE_

#define SCORE_SCALE 100000
#ifndef INF
#define INF 999999999
#endif 

struct SeedMatchType  {
  std::string seed;   // the sequence of the seed
  int pos;            // the position in the query profile where the seed is matched to
};

enum TranState {
  MM = 0, MI = 1, MD = 2, IM = 3, II = 4, DM = 5, DD = 6
};

class HMMProfile  {
 public:
  HMMProfile();
  HMMProfile(std::string &file_name);
  HMMProfile(std::list<std::string> &hmm_content);
  ~HMMProfile();
  inline int MEmitScore(int i, char a);
  inline int IEmitScore(int i, char a);
  inline int TranScore(int i, enum TranState j);
  inline int GetProfileLen(void);
  inline std::string GetName(void);
  inline char GetAAbyIndex(int i);
  std::string GetConsensus(void);
  void GenHighScoreSeeds(
      double cutoff, double seed_score_scale, ReachableReads &link_obj, 
      std::map<int, std::list<SeedMatchType> >& candidate_seeds
  );
  double CalGumblePvalue(double score);
  double CalGumbleScore(double p_value);
 protected:
  bool is_name_set_; 
  bool is_len_set_; 
  bool is_viterbi_stat_set_; 
  bool is_alphabet_set_;
  
  int profile_len_;                             // length of the profile (num. of columns)
  double viterbi_mu_;                           // (mu) for viterbi score distribution
  double viterbi_lambda_;                       // (lambda) for viterbi score distribution
  
  std::vector<char> alphabet_;                  // alphabet supported by the HMM model
  std::vector<std::vector<int> > transition_;   // transition scores from state to state
  std::vector<std::vector<int> > match_emit_;   // emission scores for match states
  std::vector<std::vector<int> > insert_emit_;  // emission scores for insertion states
  
  std::string desc_;                            // description of the profile 
  int rand_score_;                              // background score of emitting random AA
  
  // load the HMM profile
  void LoadHMM(std::list<std::string> &hmm_content);   
  // reading in corresponding fields
  bool SetName(std::string &line);
  bool SetLength(std::string &line);
  bool SetViterbiStat(std::string &line);
  bool SetScore(
      std::string &m_emit_line, 
      std::string &i_emit_line, 
      std::string &tran_line
  );
  // checking alphabet and transition settings of the HMM profile
  void CheckAlphabet(std::string &line);
  void CheckTransition(std::string &line);
  void SetRandScore(void);
  double CalMatchMerScore(int i, std::string &mer);
  
};

inline char HMMProfile::GetAAbyIndex(int i)  {
  assert(i >= 0 && i <= 19);
  switch (i)  {
    case 0: return 'A'; break;
    case 1: return 'C'; break;
    case 2: return 'D'; break;
    case 3: return 'E'; break;
    case 4: return 'F'; break;
    case 5: return 'G'; break;
    case 6: return 'H'; break;
    case 7: return 'I'; break;
    case 8: return 'K'; break;
    case 9: return 'L'; break;
    case 10: return 'M'; break;
    case 11: return 'N'; break;
    case 12: return 'P'; break;
    case 13: return 'Q'; break;
    case 14: return 'R'; break;
    case 15: return 'S'; break;
    case 16: return 'T'; break;
    case 17: return 'V'; break;
    case 18: return 'W'; break;
    case 19: return 'Y'; break;
    default: return 'X'; break;
  }
  return 'X';
}

inline int HMMProfile::MEmitScore(int i, char a) {
  // note that we might need to access the special case
  // which is encoded in the first row
  assert(i >= 0 && i <= profile_len_);  
  int j;
  switch (a)  {
    case 'A': j = 0; break;
    case 'C': j = 1; break;
    case 'D': j = 2; break;
    case 'E': j = 3; break;
    case 'F': j = 4; break;
    case 'G': j = 5; break;
    case 'H': j = 6; break;
    case 'I': j = 7; break;
    case 'K': j = 8; break;
    case 'L': j = 9; break;
    case 'M': j = 10; break;
    case 'N': j = 11; break;
    case 'P': j = 12; break;
    case 'Q': j = 13; break;
    case 'R': j = 14; break;
    case 'S': j = 15; break;
    case 'T': j = 16; break;
    case 'V': j = 17; break;
    case 'W': j = 18; break;
    case 'Y': j = 19; break;
    default: 
      //std::cout << "HMMProfile::MEmitScore: Unrecognized amino acid: " << a << std::endl;
      // output random background probability
      return rand_score_;
  }
  return match_emit_[i][j];
}

inline int HMMProfile::IEmitScore(int i, char a) {
  // note that we might need to access the special case
  // which is encoded in the first row
  assert(i >= 0 && i <= profile_len_);
  int j;
  switch (a)  {
    case 'A': j = 0; break;
    case 'C': j = 1; break;
    case 'D': j = 2; break;
    case 'E': j = 3; break;
    case 'F': j = 4; break;
    case 'G': j = 5; break;
    case 'H': j = 6; break;
    case 'I': j = 7; break;
    case 'K': j = 8; break;
    case 'L': j = 9; break;
    case 'M': j = 10; break;
    case 'N': j = 11; break;
    case 'P': j = 12; break;
    case 'Q': j = 13; break;
    case 'R': j = 14; break;
    case 'S': j = 15; break;
    case 'T': j = 16; break;
    case 'V': j = 17; break;
    case 'W': j = 18; break;
    case 'Y': j = 19; break;
    default: 
      //std::cout << "HMMProfile::IEmitScore: Unrecognized amino acid: " << a << std::endl;
      // output random background probability
      return rand_score_;
  }
  return insert_emit_[i][j];
}

inline int HMMProfile::TranScore(int i, enum TranState j) {
  // note that we might need to access the special case
  // which is encoded in the first row
  assert(i >= 0 && i <= profile_len_);
  return transition_[i][j];
}

inline int HMMProfile::GetProfileLen(void)  {
  if(!is_len_set_)  {
    std::cout << "HMMProfile::GetProfileLen: Profile length not set." << std::endl; exit(1);
  }
  return profile_len_;
}

inline std::string HMMProfile::GetName(void)  {
  if(!is_name_set_)  {
    std::cout << "HMMProfile::GetName: Profile name not set." << std::endl; exit(1);
  }
  return desc_;
}

#endif
