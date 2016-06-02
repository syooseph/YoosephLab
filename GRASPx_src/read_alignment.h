#include "reachable_reads.h"
#include "sequence_build.h"
#include "scoring_function.h"
#include "seq_align.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <deque>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <cmath>

#ifndef _READ_ALIGNMENT_
#define _READ_ALIGNMENT_

class ReadAlignment {
 public:
  ReadAlignment(void);
  ~ReadAlignment(void);
  void ComputeAlignmentScore(
      std::string& query, SequenceBuild& seq_obj,
      ScoringFunction<int>& score_obj, 
      int band_size, std::list<ReadType>& candidate_reads
  );
  void SortReadsOnEvalue(std::list<ReadType>& candidate_reads);
 protected:
  std::string GetQueryFragment(
      std::string& query, int begin, 
      int end, int band_size
  );
  
};

#endif
