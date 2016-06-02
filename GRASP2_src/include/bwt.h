#ifndef _BWT_H_
#define _BWT_H_

#include <iostream>
#include <string>
#include <list>
#include <queue>
#include <stack>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <tuple>
#include <unordered_map>

#include "divsufsort64.h"
#include "bio_alphabet.h"
#include "concatenator.h"

// the gap between two recorded occurrences in the FM-index
#define BWT_FM_GAP 64

/* Definition from 64-bit libdivsufsort */
// typedef uint8_t sauchar_t;
typedef sauchar_t BWTCHAR;
// typedef int64_t saidx_t;
typedef saidx64_t BWTIDX;
// typedef int32_t saint_t
typedef saint_t BWTINT;

class BWT {
 public:
  explicit BWT() : 
    build_success_(false),
    bwt_(NULL),
    position_(NULL),
    acc_freq_(NULL)  
  {}
  ~BWT()  {}
  void Purge(void)  {
    if(build_success_)  {
      delete [] bwt_;
      delete [] position_;
      delete [] acc_freq_;
      for(auto it = fm_index_.begin(); it != fm_index_.end(); ++ it)
        delete [] it->second;
    }
    return;
  }
  void PrintBWT(void) {std::cout << bwt_ << std::endl;}
  
  // return the size
  inline BWTIDX GetSize() {return size_;}
  
  // constructing the BWT
  void Construct(const BioAlphabet &alphabet, const char *text);
  // constructing the BWT without read ID anf position infomation
  void ConstructNoPos(const BioAlphabet &alphabet, const char *text);
  // update the current "range" while the preceeding character is "c"
  std::pair<BWTIDX, BWTIDX> UpdateRange(const char c, const std::pair<BWTIDX, BWTIDX> &range);
  // counting the occurrence of character "c" at positions in the BWT that are less than "pos" 
  inline BWTIDX CountOccurrence(const char c, const BWTIDX pos) {
    BWTIDX occ = fm_index_[c][(BWTIDX) pos / BWT_FM_GAP];
    for(BWTIDX i = pos - (pos % BWT_FM_GAP); i < pos; ++ i) occ += bwt_[i] == c ? 1 : 0;
    return occ;
  }
  // counting the occurrence of chars that are lexicographically less than "c"
  // at positions in the BWT that are less than "pos"
  BWTIDX CountLexicoLess(const char c, const BWTIDX pos);
  
  friend class BWTSearch;

 protected:
  bool build_success_;
  BioAlphabet alphabet_;  // the alphebet in use
  BWTCHAR *bwt_;          // the BWT string
  BWTIDX *position_;      // the positions of the substrings in "concat_seq"
  // the FM-index
  std::unordered_map<BWTCHAR, BWTIDX*> fm_index_;
  // the frequency of the characters     
  std::unordered_map<BWTCHAR, BWTIDX> char_freq_;
  // the accumulated frequency, which is equivalent to index in the sorted column
  BWTIDX *acc_freq_;   
  BWTIDX size_;           // the size of the text 
  
  // counting the frequency of the characters in the alphabet, excluding "delim"
  // sorts out the BWT, and also builds the FM-index
  void BuildIndex(const BWTCHAR* text);
  void BuildIndexNoPos(const BWTCHAR* text);
  
};

#endif
