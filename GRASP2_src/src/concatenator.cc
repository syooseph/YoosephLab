#include "../include/concatenator.h"

using namespace std;

Concatenator::Concatenator(char** const seq, const int n, std::string &concat_seq)  {
  if(n <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }  
  // concatenate the sequences
  concat_seq = "";
  concat_seq += DELIM;
  for(int i = 0; i < n; ++ i) {
    concat_seq += seq[i];
    concat_seq += DELIM;
  }
  return;
}

Concatenator::Concatenator(std::list<std::string> &seq, std::string &concat_seq)  {
  if(seq.size() <= 0)  {
    cerr << "Concatenator::FatalError: Empty input sequences; abort." << endl;
    exit(1);
  }  
  concat_seq = "";
  concat_seq += DELIM;
  // concatenate the sequences
  for(auto it = seq.begin(); it != seq.end(); ++ it) {
    concat_seq += *it;
    concat_seq += DELIM;
  }
  return;
}
