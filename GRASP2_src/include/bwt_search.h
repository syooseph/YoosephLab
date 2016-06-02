#ifndef _BWT_SEARCH_H_
#define _BWT_SEARCH_H_

#include <iostream>
#include <string>
#include <list>
#include <queue>
#include <stack>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <unordered_map>

#include "bwt.h"

static int cost = 10;
static int m_cost = 3;
static int g_cost = 4;
static int max_queue_size = 100;

struct ExtInfo  {
  std::pair<BWTIDX, BWTIDX> range;
  BWTINT q_pos;
  BWTINT cost;
  BWTINT cost_bound;
#ifdef DEBUG
  std::list<BWTIDX> history;
#endif
};

class ExtInfoComp {
 public:
  ExtInfoComp() {}
  
  bool operator() (const ExtInfo &lhs, const ExtInfo &rhs) const
  {
    return (lhs.cost >= rhs.cost);
  }
};

// defined as a set of forward intervals, reverse intervals, and lengths of these intervals
class IvSet {
 public:
  std::vector<std::pair<BWTIDX, BWTIDX> > ivA_, ivB_;
  std::vector<int> len_;
  
  explicit IvSet(void) {}
  void Clear() {ivA_.clear(); ivB_.clear(); len_.clear();}
  inline void PushA(const BWTIDX first, const BWTIDX second) {
    ivA_.push_back(std::make_pair(first, second));
  }
  inline void PushB(const BWTIDX first, const BWTIDX second) {
    ivB_.push_back(std::make_pair(first, second));
  }
  inline void PushLen(const int len) {
    len_.push_back(len);
  }
  bool Check(void)  {
    if(ivA_.size() != ivB_.size() || ivA_.size() != len_.size()) return false;
    for(int i = 0; i < (int) len_.size() - 1; ++ i) {
      if(len_[i] < len_[i + 1]) return false;
    }
    return true;
  }
  inline int GetSize(void)  {return len_.size();}
  void Reverse(void) {
    int n = len_.size();
    int tmp_int; std::pair<BWTIDX, BWTIDX> tmp_pair;
    for(int i = 0; i < (int) (n / 2); ++ i) {
      tmp_int = len_[i]; len_[i] = len_[n - i - 1]; len_[n - i - 1] = tmp_int;
      tmp_pair = ivA_[i]; ivA_[i] = ivA_[n - i - 1]; ivA_[n - i - 1] = tmp_pair;
      tmp_pair = ivB_[i]; ivB_[i] = ivB_[n - i - 1]; ivB_[n - i - 1] = tmp_pair;
    }
    return;
  }
};

class IvInfo  {
 public:
  BWT *bwtF_, *bwtR_;
  IvSet intervals_;
  
  explicit IvInfo(BWT *fw_bwt, BWT *re_bwt) {
    assert(fw_bwt->GetSize() == re_bwt->GetSize());
    bwtF_ = fw_bwt; bwtR_ = re_bwt;
    return;
  }
  // reset the intervals in the object while keeping the BWTs
  inline void Reset() {intervals_.Clear();}
};

struct AlignType {
  BWTIDX begin, end;
  BWTINT q_pos, cost;
};

class BWTSearch {
 public:
  explicit BWTSearch() {}
  ~BWTSearch() {}
  // search "str" against the database for exact match
  void SearchExact(BWT &bwt, const char *str, AlignType &pos);
  // a wrapper function for SearchInExact from both directions
  void Search(BWT &bwt, BWT &rev_bwt, 
      const char *str, AlignType &pos
  );
  
  // given a read of interest and minimum overlap, find all intervals (in both fw and re FM-indexes) 
  // corresponding to the prefix of the reads that perfectly overlap with the given read
  void SearchBeginIntervals(const char* seq, const int min_len, IvInfo &search_info);
  // given a set of intervals identified by SearchBeginIntervals,
  // find a set of positions that correspond to the irreducible positions in the reverse BWT
  void FindIrreducible(
      IvInfo &search_info, std::vector<BWTIDX> &ir_positions, std::vector<int> &ir_overlap
  );
  // check if the read is contained (a PROPER substring) by another read
  // if yes return true or otherwise return false 
  // also output the position in the BWT using argument "pos"
  bool IsContainedRead(const char* seq, BWT &bwt, AlignType &pos);
 protected:
  // the function calucates the lower bound and stores in the array "bound"
  // returns the length of the longest stretch of pefectly matched substring
  void CalLowerBound(BWT &rev_bwt, const char *str, int *bound);
  // search "str" against the database for in-exact mathc with up to n mismatch/gaps
  // cost is the total cost, m_cost is mismatch cost, and g_cost is gap cost
  void SearchInExact(
      BWT &bwt, BWT &rev_bwt, const int *bound,
      const char *str, AlignType &pos
  );
  
  // each phase in the searchInExact extension
  // progressively extend unless mismatch is found
  void ExtendInExactFast(
      BWT &bwt, const char *str, 
      std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
      const int *bound, int lower_cost
  );
  void Enqueue(
      std::priority_queue<ExtInfo, std::vector<ExtInfo>, ExtInfoComp> &candidate, 
      ExtInfo &phase_info
  );  
/*
*  // !!!THIS IS AN OBSOLETE FUNCTION!!!
*  // each phase in the searchInExact extension
*  void ExtendInExact(
*      BWT &bwt, const char *str, std::list<BWTIDX> &pos, 
*      std::stack<ExtInfo> &candidate, const int *bound,
*      const int cost, const int m_cost, const int g_cost    
*  );
*/
};

#endif
