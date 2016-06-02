#ifndef _MINIMIZER_SORT_H_
#define _MINIMIZER_SORT_H_

#include <iostream>
#include <unordered_map>
#include <string>
#include <cstring>
#include <list>

#include "kmer_unitcoder.h"

class MinimizerSort {
 public:
  explicit MinimizerSort() {}
  ~MinimizerSort() {}
  
  void SortSeqs(
      KmerUnitcoder &encoder, const int num_seqs, 
      char **header, char **seq
  );
  
 protected:
};

#endif
