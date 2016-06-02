#ifndef _ALIGN_BATCH_
#define _ALIGN_BATCH_

#include "scoring_prot.h"

#include <omp.h>
#include <iostream>
#include <cstring>
#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_map>


struct FullAlnType  {
  // the original query and the target sequence
  std::pair<std::string, std::string> sequences;
  // the alignment of the query and the target sequences
  std::pair<std::string, std::string> alignment;
  // the intervals that have been aligned in the query and the target sequences, respectively
  std::pair<int, int> q_interval, t_interval;
  // the alignment score
  int score;
};

class AlignBatch  {
 public:
  explicit AlignBatch(void) {}
  ~AlignBatch(void) {}
  
  // globally aligning the pivot sequence "sp" to "n" sequences in "sgroup"
  // band-size is set to "b" with gap open "g", extension "e", and scoring function "f"
  // the alignment scores are stored in "r" 
  void AlignGlobal(
      char *sp, const int n, char **sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  void AlignGlobal(
      std::string &sp, const int n, std::vector<std::string> &sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  // compute pairwise alignment and store the score in r
  // return the resulting alignment score
  int AlignGlobalPairwise(
      const std::string &seq1, const std::string &seq2, 
      const int band, ScoringProt &f
  );
  // globally aligning two sequences "s" and "t" using the edit distance scoring function
  // return the edit distance
  int AlignGlobalSingleEditDist(const std::string &s, const std::string &t, int band);
  
  int AlignGlobalSingleEditDist(
      const std::string &s, const std::string &t, 
      ScoringProt &f, const int gap_cost, int band
  );
  
  // align pairwise single sequence completely, use larger memory and longer time
  // but will return the exact aligned regions in the query and the target
  // and will also generate the actual alignment
  // returns the alignment score
  int AlignGlobalSingle(
      const std::string &seq1, const std::string &seq2, 
      ScoringProt &f, std::pair<std::string, std::string> &alignment
  );
  
  // locally aligning the pivot sequence "sp" to "n" sequences in "sgroup"
  // band-size is set to "b" with gap open "g", extension "e", and scoring function "f"
  // the alignment scores are stored in "r" 
  void AlignLocal(
      char *sp, const int n, char **sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  void AlignLocal(
      std::string &sp, const int n, std::vector<std::string> &sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  // compute pairwise alignment and store the score in r
  // return the resulting alignment score
  int AlignLocalPairwise(
      const std::string &seq1, const std::string &seq2, 
      const int band, ScoringProt &f
  );
  
  // align pairwise single sequence completely, use larger memory and longer time
  // but will return the exact aligned regions in the query and the target
  // and will also generate the actual alignment
  // returns the alignment score
  int AlignLocalSingle(
      const std::string &seq1, const std::string &seq2, 
      ScoringProt &f, std::pair<std::string, std::string> &alignment,
      std::pair<int, int> &q_interval, std::pair<int, int> &t_interval
  );
  
  // using multi-threading to parallelize the global alignment
  void MultiAlignGlobal(
      const int threads, char *sp, const int n, char **sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  void MultiAlignGlobal(
      const int threads, std::string &sp, const int n, std::vector<std::string> &sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  // using multi-threading to parallelize the global alignment
  // aligning n pairs of sequences, each pair of sequences are aligned
  void MultiAlignGlobalPairwise(
      const int threads, 
      std::vector<std::string> &sgroup1, std::vector<std::string> &sgroup2,
      const int band, ScoringProt &f,std::vector<int> &r
  );
  
  // using multi-threading to parallelize the local alignment
  // aligning one sequence against a collection of sequences
  void MultiAlignLocal(
      const int threads, char *sp, const int n, char **sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  void MultiAlignLocal(
      const int threads, std::string &sp, const int n, std::vector<std::string> &sgroup, 
      const int band, ScoringProt &f, std::vector<int> &r
  );
  // perform the alignment and trace-back (printing information available using this function)
  void MultiAlignLocalFull(
      const int threads, std::string &sp, const int n, std::vector<std::string> &sgroup, 
      ScoringProt &f, std::vector<FullAlnType> &r
  );
  // using multi-threading to parallelize the local alignment
  // aligning n pairs of sequences, each pair of sequences are aligned
  void MultiAlignLocalPairwise(
      const int threads, 
      std::vector<std::string> &sgroup1, std::vector<std::string> &sgroup2,
      const int band, ScoringProt &f, std::vector<int> &r
  );
  
  // compute the kmer similarity between the two sequences
  // n is the mer length, seq1 and seq2 are the two sequences
  int KmerSimilarity(const int n, std::string &seq1, std::string &seq2);
  // compute the kmer similarity between the two sequences
  // n is the mer length, seq1 and seq2 are the two sequences
  int WeightedKmerSimilarity(const int n, ScoringProt &f, std::string &seq1, std::string &seq2);
  
  int Max3(int a, int b, int c);
  int Min3(int a, int b, int c);
};

#endif
