#include <boost/algorithm/string.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>
#include <cstdlib>

#define SCORE_SCALE 100000

class HMMProfile  {
 public:
  HMMProfile();
  HMMProfile(std::string &file_name);
  ~HMMProfile();
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
  
  // load the HMM profile
  void LoadHMM(std::list<std::string> &hmm_content);   
  // reads in contents in the NAME field
  inline bool SetName(std::string &line) {
    std::vector<std::string> decom;
    boost::split(decom, line, boost::is_any_of("\t "));
    if(decom[0] == "NAME" && decom.size() == 2)  {
      desc_ = decom[1];
      is_name_set_ = true;
      return true;
    }
    return false;
  }
};
