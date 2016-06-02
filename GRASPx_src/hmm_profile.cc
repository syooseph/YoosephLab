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
}

HMMProfile::~HMMProfile() {
  ;
}

void HMMProfile::LoadHMM(std::string &file_name)  {
  // parse the header information
  auto it = hmm_content.begin();
  if(it->substr(0, 6) != string("HMMER3"))  {
    cout << "HMMProfile::LoadHMM: Error: Unrecognized HMM profile file format." << endl;
    exit(1);
  }
  ++ it;  // go to the next line
  
  for(it; it != hmm_content.end(); ++ it) {
    // need NAME, LENG, HMM, STATS LOCAL VITERBI
    if(it->substr(0, 4) == string("NAME"))  {
      is_name_set_ = SetName(*it);;
    } else if(it->substr(0, 4) == string("LENG")) {
      is_len_set_ = SetLength(*it);;
    } else if(it->substr(0, 19) == string("STATS LOCAL VITERBI")) {
      is_vitrbi_stat_set_ = SetViterbiStat(*it);;
    } else if(it->substr(0, 3) == string("HMM")) {
      is_alphabet_set_ = SetAlphbet(*it);;
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
    ++ it;  // skip the transition state label, assuming the following order
            // "m->m     m->i     m->d     i->m     i->i     d->m     d->d"
  }
  // now the iterator should be set on the first line of the scores
  // need to check if it corresponds the 
  while(it != hmm_content.end() && *it != string("//")) {
    // the scores come in 3 together
    string m_emit_line, i_emit_line, tran_line;
    if(it != hmm_content.end()) {m_emit_line = *it; ++ it};
    if(it != hmm_content.end()) {i_emit_line = *it; ++ it};
    if(it != hmm_content.end()) {tran_line = *it; ++ it};
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
    if(transition_.size() != profile_len_ ||
       match_emit_.size() != profile_len_ ||
       insert_emit_.size() != profile_len_
    )  {
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
