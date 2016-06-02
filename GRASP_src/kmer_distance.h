#include "index_sample.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <string>
#include <list>
#include <queue>

#ifndef _KMER_DISTANCE
#define _KMER_DISTANCE

struct KmerIterType {
  std::string kmer;
  int iter;
};

class KmerDistance {
 public:
  KmerDistance(void);
  KmerDistance(unsigned int in_kmer_size, enum Alphabet in_alphabet);
  ~KmerDistance(void);  
  void KmerAdjacency(std::string file_in); // define adjacency between k-mers 
  void GetKmerClosure(
    const double& conv_rate, 
    const std::unordered_map<std::string, bool>& seed_kmers, 
    std::unordered_map<std::string, bool>& kmer_closure
  );  // define kmer_closure
 private:
  // shared variables
  std::unordered_map<std::string, bool> boundary_kmers_;
  std::unordered_map<std::string, std::unordered_map<std::string, bool> > kmer_adj_list_;
  std::queue<KmerIterType> candidate_kmer_;
  unsigned int kmer_size_;
  enum Alphabet alphabet_;
  // functions
  void UpdateAdjacency(std::string seq);
};

#endif
