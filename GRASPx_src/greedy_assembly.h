#include "reachable_reads.h"
#include "sequence_build.h"

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

#ifndef _GREEDY_ASSEMBLY_
#define _GREEDY_ASSEMBLY_

typedef std::list<ReadType>::iterator ReadListIterType;

class GreedyAssembly  {
 public:
  GreedyAssembly(void);
  ~GreedyAssembly(void);
  void BuildRIDLink(
      std::list<ReadType>& candidate_reads, 
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
  );
  void GreedyExtend(
      ReachableReads& index_obj,
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link,
      ReadType seed_read, double evalue_cutoff,
      std::deque<ReadConnectType>& contig_path
  );
  std::string SpellContigSequence(
      SequenceBuild& seq_obj, std::deque<ReadConnectType>& contig_path
  );
 protected:
  ReadConnectType GetBestConnectReadFW(
      ReachableReads& index_obj, ReadConnectType& ext_read, 
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
  );
  ReadConnectType GetBestConnectReadRE(
      ReachableReads& index_obj, ReadConnectType& ext_read, 
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link
  );
  void FWAssemble(
      ReachableReads& index_obj, 
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link, 
      double evalue_cutoff, std::deque<ReadConnectType>& contig_path
  );
  void REAssemble(
      ReachableReads& index_obj, 
      std::unordered_map<RIDType, std::list<ReadListIterType> >& rid_link, 
      double evalue_cutoff, std::deque<ReadConnectType>& contig_path
  );
  bool IsPreRegionOverlap(int r1_begin, int r1_end, int r2_begin, int r2_end);
  bool IsAftRegionOverlap(int r1_begin, int r1_end, int r2_begin, int r2_end);
};

#endif
