#ifndef _ITERATIVE_ASSEMBLY_H_
#define _ITERATIVE_ASSEMBLY_H_

#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "debruijn_graph.h"
#include "bwt.h"
#include "bwt_search.h"
#include "kmer_unitcoder.h"
#include "frequency_table.h"


struct UnitigType {
  std::string seq;
};

class IterativeAssembly {
 public:
  // num_seqs: the number of raw reads
  // seqs: the actual two-dimensional array that actually holds the sequences
  // kmer_config: a set of k-mers to be tried out for the iterative assembly
  // coverage_config: a set of coverages to be tried out for the iterative assembly
  explicit IterativeAssembly(
      BioAlphabet &alphabet,
      int num_seqs, char **seqs, 
      std::vector<int> &kmer_config, std::vector<int> &coverage_config
  ) : alphabet_(alphabet),
      num_seqs_(num_seqs), 
      seqs_(seqs), 
      kmer_config_(kmer_config), 
      coverage_config_(coverage_config) 
  {
    if(!CheckConfig())  {
      std::cerr << "Impropertly set k-mer sizes or coverages... Abort." << std::endl;
      exit(1);  
    }
  }
  ~IterativeAssembly()  {}
  // Driver function for the iterative assembly
  // unitigs: the output of the iterative assembly, a set of unitigs with reads mapped on them
  // read_mark: marking whether a read gets mapped to some of the unitigs
  void Assemble(std::vector<UnitigType> &unitigs, std::vector<bool> &read_mark);
 private:
  BioAlphabet alphabet_;
  int num_seqs_;
  char **seqs_;
  std::vector<int> kmer_config_; 
  std::vector<int> coverage_config_;
  // checking whether the configuration is set correctly 
  bool CheckConfig();
  // counting frequency filtering table
  void ConstructFreqTable(
      const int mer_len, 
      std::vector<bool> &read_mark, 
      FrequencyTable<KmerUnitType> &unit_freq
  );
  // subtracting the k-mer frequencies for the recruited reads
  void UpdateFreqTable(
      const int mer_len, 
      std::vector<bool> &read_mark, 
      FrequencyTable<KmerUnitType> &unit_freq
  );
};

#endif
