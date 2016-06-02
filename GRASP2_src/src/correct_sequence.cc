#include "../include/correct_sequence.h"

using namespace std;

CorrectSequence::CorrectSequence()  {
  return;
}

CorrectSequence::~CorrectSequence() {
  return;
}

int CorrectSequence::CorrectError(
    char **seq, int n, BioAlphabet &alphabet, 
    int mer_len, KmerGraph &kmer_freq, std::set<int> &uncorrectable,
    int max_iter
) {
  //cout << "maximum iteration: " << max_iter << endl;
  // boundary case initialization
  CInfoType init;   // the recent string and corrections are empty 
  init.cpos = mer_len - 2;
  init.score = 0.0;
  KmerEncoder encoder(alphabet, mer_len);
  // the IDs of the sequences that we failed to correct
  list<int> failed_corrections;
  int num_success = 0;
  for(int i = 0; i < n; ++ i) {
    //cout << seq[i] << endl;
    // for sequences that are too short
    if(strlen(seq[i]) < mer_len)  {
      failed_corrections.push_back(i);
      continue;
    }
    // initialize the solution pool
    priority_queue<CInfoType> solution_pool;
    init.recent = 'X' + string(seq[i]).substr(0, mer_len - 1);  // 'X' is a leading foo character
    solution_pool.push(init);
    // implementing Heng Li's BFC algorithm
    int num_iter = 0;
    CInfoType best_solution;
    best_solution.cpos = -1;
    best_solution.score = (float) -MAX_SCORE;
    while(!solution_pool.empty() && num_iter < max_iter) {
      // get the first element in the queue
      CInfoType current = solution_pool.top();
      solution_pool.pop();
      //cout << num_iter << " " << current.cpos << endl;
      if(current.cpos >= strlen(seq[i]) - 1)  {
        //cout << "reaching the end" << endl;
        // a solution has been computed
        if(current.score > best_solution.score) {
          best_solution = current;
          break;
        }
      } else  {
        string recent_sx = current.recent.substr(1, mer_len - 1);   
        string recent_nw = recent_sx + seq[i][current.cpos + 1];  
        //cout << "Current kmer:  " << recent_nw << endl;  
        // a solution needs to be further explored  
        // note that one should remove low-coverage mers using KmerFrequency::RemoveLowFreqKmer
        if(kmer_freq.GetKmerFreq(encoder.Encode(recent_nw.c_str())) > 0)  { 
          //cout << "---No need correction" << endl;
          // this is a high frequency k-mer, do not need to correct it
          ++ current.cpos;
          current.recent = recent_nw;
          solution_pool.push(current);
        } else  {
          //cout << "---Need correction" << endl;
          for(auto it_a = alphabet.char_map_.begin(); it_a != alphabet.char_map_.end(); ++ it_a) {
            recent_nw = recent_sx + it_a->first;
            int mer_freq = kmer_freq.GetKmerFreq(encoder.Encode(recent_nw.c_str()));
            if(mer_freq > 0)  {
              //cout << "---Corrected" << endl;
              CInfoType c_next = current;
              ++ c_next.cpos;
              c_next.score -= 1.0 + (1.0 / mer_freq);
              c_next.recent = recent_nw;
              c_next.corrections.push_back(pair<int, char>(c_next.cpos, it_a->first));
              solution_pool.push(c_next);
            }
          }
        }
      }
      ++ num_iter;
    }
    // correct the sequence
    if(best_solution.cpos >= strlen(seq[i]) - 1 && best_solution.score > (float) -MAX_SCORE)  {   
      for(
          auto it_c = best_solution.corrections.begin(); 
          it_c != best_solution.corrections.end(); ++ it_c
      ) {
        // make correction to the sequence
        seq[i][it_c->first] = it_c->second;
      }
      ++ num_success;
    } else  {
      uncorrectable.insert(i);
    }
  }
  return num_success;
}
