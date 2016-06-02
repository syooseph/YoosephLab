#include "gsa.h"
#include "sequence.h"

#include "timer.h"
#include <boost/filesystem.hpp>
#include <cstring>
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <map>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <cstdlib>
#include <set>

using namespace std;

string rev_comp(const string& seq) {
  string rc_seq = "";
  for(auto it = seq.rbegin(); it != seq.rend(); ++ it) {
    char c = *it;
    if(c == 'A')  {
      rc_seq += 'U';
    } else if(c == 'U') {
      rc_seq += 'A';
    } else if(c == 'C') {
      rc_seq += 'G';
    } else if(c == 'G') {
      rc_seq += 'C';
    }
  }
  return rc_seq;
}

int main()  {
  string seq = "ACUACGUCGUCUGAUCGUACGUAGCGCGUACGAUCGUACGACGUAAUUAGCGCGCUAUAGGCAUGCAUCGAGCGCAUCGUAAGGCUCGACAGUCACUACGUCGUCUGAUCGUACGUAGCGCGUACGAUCGUACGACGUAAUUAGCGCGCUAUAGGCAUGCAUCGAGCGCAUCGUAAGGCUCGACAGUC";
  string rc_seq = rev_comp(seq);
  char **sample_seqs = new char*[2];
  sample_seqs[0] = new char[seq.length() + 1];
  strcpy(sample_seqs[0], seq.c_str());
  sample_seqs[1] = new char[rc_seq.length() + 1];
  strcpy(sample_seqs[1], rc_seq.c_str());
  GSA* suffix_array = new GSA(sample_seqs, 2, true);
  
  int begin = 0, end = 0;
  int index;
  int len_cutoff = 3;
  list<pair<int, int> > palindroms;
  for(index = 0; index < suffix_array->getSize() - 1; ++ index) {
    if(suffix_array->getId(index) != suffix_array->getId(index + 1))  {
      int len = suffix_array->getLcpAt(index + 1);
      if(len < len_cutoff)  {
        continue;
      }
      string left, right;
      if(suffix_array->getId(index) == 0)  {
        left = seq.substr(suffix_array->getPos(index), len);
        right = seq.substr(seq.length() - suffix_array->getPos(index + 1) - len, len);
      } else  {
        left = seq.substr(suffix_array->getPos(index + 1), len);
        right = seq.substr(seq.length() - suffix_array->getPos(index) - len, len);
      }
      cout << left << " " << right << endl;
    }
  }
  return 0;
}
