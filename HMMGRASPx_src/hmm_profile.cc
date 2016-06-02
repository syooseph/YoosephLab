#include "hmm_profile.h"

using namespace std;

HMMProfile::HMMProfile()  {
  is_name_set_ = is_len_set_ = is_viterbi_stat_set_ = is_alphabet_set_ = false;
}

HMMProfile::HMMProfile(std::string &file_name)  {
  is_name_set_ = is_len_set_ = is_viterbi_stat_set_ = is_alphabet_set_ = false;
  ifstream HMMFILE(file_name);
  if(!HMMFILE.is_open())  {
    cout << "HMMProfile::LoadHMM: Error: Cannot open HMM profile file." << endl;
    exit(1);
  }
  // reads in the HMM info
  list<string> hmm_content;
  string line;
  while(getline(HMMFILE, line)) {
    hmm_content.push_back(line);
  }
  HMMFILE.close();
  LoadHMM(hmm_content);
  SetRandScore();
}

HMMProfile::HMMProfile(std::list<std::string> &hmm_content) {
  is_name_set_ = is_len_set_ = is_viterbi_stat_set_ = is_alphabet_set_ = false;
  LoadHMM(hmm_content);
  SetRandScore();
}

HMMProfile::~HMMProfile() {
  ;
}


void HMMProfile::LoadHMM(std::list<std::string> &hmm_content)  {
  // parse the header information
  auto it = hmm_content.begin();
  if(it->substr(0, 6) != string("HMMER3"))  {
    cout << "HMMProfile::LoadHMM: Error: Unrecognized HMM profile file format." << endl;
    exit(1);
  }
  ++ it;  // go to the next line
  
  for(; it != hmm_content.end(); ++ it) {
    // need NAME, LENG, HMM, STATS LOCAL VITERBI
    if(it->substr(0, 4) == string("NAME"))  {
      is_name_set_ = SetName(*it);;
    } else if(it->substr(0, 4) == string("LENG")) {
      is_len_set_ = SetLength(*it);;
    } else if(it->substr(0, 19) == string("STATS LOCAL VITERBI")) {
      is_viterbi_stat_set_ = SetViterbiStat(*it);;
    } else if(it->substr(0, 3) == string("HMM")) {
      CheckAlphabet(*it);;
      break;
    }
  }
  // requires alphabet to be set
  if(!is_alphabet_set_)  {
    cout << "HMMProfile::LoadHMM: Error: Fail to load amino alphabet (corrputed file or HMM is built on nucleotide space)." << endl;
    exit(1);
  }
  // go for recording the scores
  if(it != hmm_content.end())  {
    ++ it;  // go to the transition header line
    // Assuming the following order
    // "m->m     m->i     m->d     i->m     i->i     d->m     d->d"
    CheckTransition(*it);
    ++ it;  // go to the scores  
  }
  // now the iterator should be set on the first line of the scores
  // need to check if it corresponds the background COMPO field
  if(it->find("COMPO") != string::npos)  {
    string m_line = *it; ++ it;
    string i_line = *it; ++ it;
    string t_line = *it; ++ it;
    SetScore(m_line, i_line, t_line);
  } else  {
    cout << "HMMProfile::LoadHMM: Error: Cannot find background composition in profile." << endl;
    exit(1);
  }
  while(it != hmm_content.end() && *it != string("//")) {
    // the scores come in 3 together
    string m_emit_line, i_emit_line, tran_line;
    if(it != hmm_content.end()) {m_emit_line = *it; ++ it;};
    if(it != hmm_content.end()) {i_emit_line = *it; ++ it;};
    if(it != hmm_content.end()) {tran_line = *it; ++ it;};
    bool set_succ;
    if(!m_emit_line.empty() && !i_emit_line.empty() && !tran_line.empty())  {
      set_succ = SetScore(m_emit_line, i_emit_line, tran_line);
    }
    if(!set_succ)  {
      cout << "HMMProfile::LoadHMM: Error: Fail to set scores. HMM profile file might be corrupted." << endl;
      exit(1);
    }
  }
  // check for number of columns, see if matches what defined by the header
  if(is_len_set_)  {
    if(transition_.size() != (unsigned int) profile_len_ + 1 ||
       match_emit_.size() != (unsigned int) profile_len_ + 1||
       insert_emit_.size() != (unsigned int) profile_len_ + 1
    )  {
      cout << transition_.size() << " " << match_emit_.size() << endl;
      cout << "HMMProfile::LoadHMM: Error: Number of scores do not match pre-defined sequence length. HMM profile file might be corrupted." << endl;
      exit(1);
    }
  } else  {
    if(transition_.size() == match_emit_.size() ||
       match_emit_.size() == insert_emit_.size()
    )  {
      profile_len_ = transition_.size();
      is_len_set_ = true;
    } else  {
      cout << "HMMProfile::LoadHMM: Error: Number of scores do not match pre-defined sequence length. HMM profile file might be corrupted." << endl;
      exit(1);
    }
  }
  return;
}

// reads in contents in the NAME field
bool HMMProfile::SetName(std::string &line) {
  std::vector<std::string> decom;
  boost::split(decom, line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << line << " " << decom[0] << "  " << decom[1] << endl;
  if(decom[0] == "NAME" && decom.size() == 2)  {
    desc_ = decom[1];
    //cout << desc_ << endl;
    is_name_set_ = true;
    return true;
  }
  return false;
}

void HMMProfile::SetRandScore(void) {
  // the first row of the emission score describes the background emission scores
  double s = 0;
  for(int i = 0; i < 20; ++ i) {
    s += match_emit_[0][i];
  }
  s /= 20;
  rand_score_ = (int) s * SCORE_SCALE;
  return;
}

// reads in the contents in the LENG field
bool HMMProfile::SetLength(std::string &line) {
  std::vector<std::string> decom;
  boost::split(decom, line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << line << " " << decom[0] << "  " << decom[1] << endl;
  if(decom[0] == "LENG" && decom.size() == 2)  {
    profile_len_ = std::stoi(decom[1]);
    //cout << profile_len_ << endl;
    is_len_set_ = true;
    return true;
  }
  return false;
}

bool HMMProfile::SetViterbiStat(std::string &line) {
  std::vector<std::string> decom;
  boost::split(decom, line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << line << " " << decom[0] << "  " << decom[1] << endl;
  if(decom[0] == "STATS" && decom[1] == "LOCAL" && 
     decom[2] == "VITERBI" && decom.size() == 5
  )  {
    viterbi_mu_ = std::stof(decom[3]);
    viterbi_lambda_ = std::stof(decom[4]);
    //cout << viterbi_mu_ << "  " << viterbi_lambda_ << endl;
    is_viterbi_stat_set_ = true;
    return true;
  }
  return false;
}

bool HMMProfile::SetScore(
    std::string &m_emit_line, 
    std::string &i_emit_line, 
    std::string &tran_line
) {
  // check for current scoring matrix consistency
  if(!(match_emit_.size() == insert_emit_.size() && 
       insert_emit_.size() == transition_.size())
  )  {
    cout << "HMMProfile::SetScore: Error: Incompatible sizes of scoring matrices." << endl;
    exit(1);
  }
  std::vector<std::string> decom_m, decom_i, decom_t;
  boost::split(decom_m, m_emit_line, boost::is_any_of("\t "), boost::token_compress_on);
  boost::split(decom_i, i_emit_line, boost::is_any_of("\t "), boost::token_compress_on);
  boost::split(decom_t, tran_line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << m_emit_line << endl;
  //cout << i_emit_line << endl;
  //cout << tran_line << endl;
  //cout << decom_m.size() << " " << decom_i.size() << "  " << decom_t.size() << endl;
  if(decom_m.size() >= 22 && decom_i.size() >= 21 && decom_t.size() >= 8)  {
    // loads in scores, assuming alphabet size is 20 (amino acids)
    vector<int> score_m(20, 0);
    for(unsigned int i = 2; i < 22; ++ i) {
      if(decom_m[i] == "*") score_m[i - 2] = INF;
      else score_m[i - 2] = (int) (std::stof(decom_m[i]) * SCORE_SCALE); 
    }
    vector<int> score_i(20, 0);
    for(unsigned int i = 1; i < 21; ++ i) {
      if(decom_i[i] == "*") score_i[i - 1] = INF;
      else score_i[i - 1] = (int) (std::stof(decom_i[i]) * SCORE_SCALE); 
    }
    vector<int> score_t(7, 0);
    for(unsigned int i = 1; i < 8; ++ i) {
      if(decom_t[i] == "*") score_t[i - 1] = INF;
      else score_t[i - 1] = (int) (std::stof(decom_t[i]) * SCORE_SCALE); 
    }
    // record the scores
    match_emit_.push_back(score_m);
    insert_emit_.push_back(score_i);
    transition_.push_back(score_t);
    // check for sizes. note that we need to consider special transitions
    if(decom_m[1] != "COMPO")  {
      int state_index = std::stoi(decom_m[1]);
      assert(match_emit_.size() == (unsigned int) state_index + 1);
      assert(insert_emit_.size() == (unsigned int) state_index + 1);
      assert(transition_.size() == (unsigned int) state_index + 1);
    }
    return true;
  } else  {
    cout << "HMMProfile::SetScore: Error: Incompatible sizes of scoring matrices. HMM profile file might be corrupted." << endl;
    exit(1);
  }
  return false;
}

void HMMProfile::CheckAlphabet(std::string &line) {
  std::vector<std::string> decom;
  boost::split(decom, line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << line << endl;
  if(decom.size() >= 21 && 
     decom[1] == "A" && decom[2] == "C" && decom[3] == "D" && decom[4] == "E" &&
     decom[5] == "F" && decom[6] == "G" && decom[7] == "H" && decom[8] == "I" &&
     decom[9] == "K" && decom[10] == "L" && decom[11] == "M" && decom[12] == "N" &&
     decom[13] == "P" && decom[14] == "Q" && decom[15] == "R" && decom[16] == "S" &&
     decom[17] == "T" && decom[18] == "V" && decom[19] == "W" && decom[20] == "Y"
   )  {
    is_alphabet_set_ = true;
  } else  {
    cout << "HMMProfile::CheckAlphabet: Error: Amino acids arranged out of order (alphabetical order of amino acids expected)." << endl;
    exit(1);
  }
  return;
}

void HMMProfile::CheckTransition(std::string &line)  {
  std::vector<std::string> decom;
  boost::split(decom, line, boost::is_any_of("\t "), boost::token_compress_on);
  //cout << line << endl;
  if(decom.size() >= 8 && 
     decom[1] == "m->m" && decom[2] == "m->i" && decom[3] == "m->d" &&
     decom[4] == "i->m" && decom[5] == "i->i" && decom[6] == "d->m" &&
     decom[7] == "d->d"
  )  {
    return;
  } else  {
    cout << "HMMProfile::CheckAlphabet: Error: Transition states arranged out of order (expected \"m->m     m->i     m->d     i->m     i->i     d->m     d->d\")." << endl;
    exit(1);
  }
  return;
}

std::string HMMProfile::GetConsensus()  {
  string consensus;
  for(int i = 1; i <= profile_len_; ++ i) {
    int best_j = 0;
    double best_score = match_emit_[i][0];
    for(int j = 1; j < 20; ++ j)  {
      if(match_emit_[i][j] < best_score)  {
        best_score = match_emit_[i][j];
        best_j = j;
      }
    }
    consensus += GetAAbyIndex(best_j);
  }
  return consensus;
}

void HMMProfile::GenHighScoreSeeds(
    double cutoff, double seed_score_scale, ReachableReads &link_obj,
    std::map<int, std::list<SeedMatchType> >& candidate_seeds
)  {
  
  assert(profile_len_ >= link_obj.seed_len_);
  DatabaseIndex db_obj_foo;
  string query = GetConsensus();
  char *query_str = new char [query.length() + 1];
  strcpy(query_str, query.c_str());
  map<int, list<SeedType> > seed_holder;
  // create a list of position where each 3-mer can be used as a high-score mer
  map<uint16_t, map<int, bool> > cand_q_pos;
  for(int k = 0; k < (int) query.length() - 3; ++ k)  {
    uint16_t q_idx = db_obj_foo.Convert3mer(&query_str[k]);
    auto it = link_obj.high_score_match_.find(q_idx);
    if(it != link_obj.high_score_match_.end())  {
      for(auto it_m = it->second.begin(); it_m != it->second.end(); ++ it_m)  {
        cand_q_pos[*it_m][k] = true;
      }
    }
  }
  // check for each presented seeds and locate their position in profile
  for(auto it_s = link_obj.seed_ext_.begin(); it_s != link_obj.seed_ext_.end(); ++ it_s)  {
    uint16_t head_index = db_obj_foo.Convert3mer(it_s->first.substr(0, 3).c_str());
    uint16_t tail_index = db_obj_foo.Convert3mer(it_s->first.substr(link_obj.seed_len_ - 3, 3).c_str());
    if(cand_q_pos.find(head_index) == cand_q_pos.end() || cand_q_pos.find(tail_index) == cand_q_pos.end())  {
      continue;
    }
    double best_match_score = INF;
    int best_position = -1;
    for(auto it_pre = cand_q_pos[head_index].begin(); it_pre != cand_q_pos[head_index].end(); ++ it_pre) {
      //cout << i << endl;
      if(cand_q_pos[tail_index].find(it_pre->first + link_obj.seed_len_ - 3) != cand_q_pos[tail_index].end())  {
        string seed_mer = it_s->first;
        double match_score = CalMatchMerScore(it_pre->first + 1, seed_mer);
        if(match_score < best_match_score)  {
          best_match_score = match_score;
          best_position = it_pre->first;
        }
      }
    }
    if(best_match_score < cutoff * SCORE_SCALE && best_position >= 0)  {
      SeedType sd;
      sd.seed_seq = it_s->first;
      sd.q_pos = best_position + 1;
      seed_holder[(int) best_match_score].push_back(sd); 
    }
  }
  // remove seed redundancies
  link_obj.RefineSeed(ceil(seed_score_scale * link_obj.seed_len_), query.length(), seed_holder);
  // format into seedpair type
  for(auto it = seed_holder.rbegin(); it != seed_holder.rend(); ++ it) {
    for(auto it_s = it->second.begin(); it_s != it->second.end(); ++ it_s) {
      if(link_obj.seed_ext_.find(it_s->seed_seq) != link_obj.seed_ext_.end())  {
        SeedMatchType seed_match;
        seed_match.seed = it_s->seed_seq;
        seed_match.pos = it_s->q_pos;
        candidate_seeds[it->first].push_back(seed_match);
      }
    }
  }
  return;
}

double HMMProfile::CalMatchMerScore(int i, std::string &mer) {
  assert(i >= 1 && i + (int) mer.length() - 1 <= profile_len_);
  double score = 0;
  for(unsigned int k = 0; k < mer.length(); ++ k) {
    // assuming transition score is 0
    // add emit score
    score += MEmitScore(i + k, mer[k]);
    // subtract background score
    score -= MEmitScore(0, mer[k]);
  }
  // note that the score returned is the exponent of the probability that the sequence
  // mer is observed at random at position i of the profile
  return score;
}

double HMMProfile::CalGumblePvalue(double score) {
  double x = -((score / SCORE_SCALE) - viterbi_mu_) / viterbi_lambda_;
  if(x < -5)  {
    return 0;
  } else if(x > 5) {
    return INF;
  } else  {
    return exp(x);
  }
}

double HMMProfile::CalGumbleScore(double p_value) {
  //cout << "Viterbi info:  " << viterbi_mu_ << " " << viterbi_lambda_ << endl;
  double x = -log(p_value);
  return viterbi_mu_ + viterbi_lambda_ * log(x);
}
