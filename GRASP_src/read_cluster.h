#include "index_sample.h"
#include "gsa.h"
#include "file.h"
#include "sequence.h"

#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <queue>
#include <list>

#ifndef _READ_CLUSTER_
#define _READ_CLUSTER_

class ReadCluster {
 public:
  ReadCluster(
      const int& n_back, const std::string& file_name, const std::string& lcp_file, 
      const std::string& mcp_file, const std::string& gsa_file
  );
  ~ReadCluster();
  void RecruitConnectedReads(void);
  void InterpretConnections(void);
  void GetClosureSequences(const std::list<std::string>& seed_kmers);
  std::string GetReadSequence(const int& rank);
 private:
  GSA *suffix_array_;
  char **sequence_;
  int num_sequences_;
  std::vector<bool> reads_resolved;
  std::unordered_map<KmerType, int> kmer_rank_;
  std::unordered_map<int, std::unordered_map<int, bool> > connected_kmers_;
  std::unordered_map<int, int> kmer_cluster_; // read_ID -> cluster_ID
  std::unordered_map<int, bool> recruited_kmer_;  // read_ID -> T/F
  int n_back_len_;
  void LoadSequence(const std::string& file_name);
  void LoadSuffixArray(
      const std::string& lcp_file, const std::string& mcp_file, const std::string& gsa_file
  );
};

#endif
