#include "../include/minimizer_sort.h"

using namespace std;

void MinimizerSort::SortSeqs(
    KmerUnitcoder &encoder, const int num_seqs, 
    char **header, char **seq
)  { 
  // make a copy of the sequences
  char **header_holder = new char* [num_seqs];
  char **seq_holder = new char* [num_seqs];
  map<KmerUnitType, list<int> > orders;
  int i, j;
  int mer_len = encoder.GetMerLen();
  
  for(i = 0; i < num_seqs; ++ i) {
    header_holder[i] = header[i];
    seq_holder[i] = seq[i];
    KmerUnitType min_en = encoder.Encode(seq[i]);
    KmerUnitType en = min_en;
    for(j = mer_len; j < strlen(seq[i]); ++ j) {
      en = encoder.RightExt(en, seq[i][j]);
      if(en < min_en) min_en = en;
    }
    orders[min_en].push_back(i);
    //cout << "The " << i << "th sequence done..." << endl;
  }
  // copy the sequences
  int nn = 0;
  for(auto it = orders.begin(); it != orders.end(); ++ it) {
    for(auto it_l = it->second.begin(); it_l != it->second.end(); ++ it_l) {
      header[nn] = header_holder[*it_l];
      seq[nn] = seq_holder[*it_l];
      ++ nn;
    }
  }
  // release the memory
  delete [] header_holder;
  delete [] seq_holder;
  return;
}
