#ifndef _CORRECT_SEQUENCE_H_
#define _CORRECT_SEQUENCE_H_

#include <iostream>
#include <string>
#include <cstring>
#include <tuple>
#include <queue>
#include <set>

#include "kmer_encoder.h"
#include "bio_alphabet.h"
#include "kmer_graph.h"

#ifndef MAX_SCORE
#define MAX_SCORE 1000000
#endif

// the first integer is the index of the sequence
// the second char is the character to replace the original sequence
typedef std::pair<int, char> CorrectionType;

class CInfoType  {
 public:
  int cpos;   // the current position
  float score;   // the cost so-far     
  std::string recent;   // the most recent k-long string
  std::list<CorrectionType> corrections;  // solution for the current correction
  // overloading the operators for priority queue operations
  friend bool operator> (const CInfoType &c1, const CInfoType &c2) {
    return (c1.score > c2.score);
  }
  friend bool operator>= (const CInfoType &c1, const CInfoType &c2) {
    return (c1.score >= c2.score);
  }
  friend bool operator< (const CInfoType &c1, const CInfoType &c2) {
    return (c1.score < c2.score);
  }
  friend bool operator<= (const CInfoType &c1, const CInfoType &c2) {
    return (c1.score <= c2.score);
  }
};

class CorrectSequence {
 public:
  CorrectSequence();
  ~CorrectSequence();
  // Correct errors for sequences in "seq"
  // returns the number of corrected sequences (expected to be equal or less than n)
  // seq: two-dimensional table holding the sequences
  // n: number of sequences in seq
  // alphabet: sequence alphabet for seq
  // mer_len: the length of the k-mer
  // kmer_freq: the frequency of all k-mers
  // max_iter:  maximum number of allowed iterations
  int CorrectError(
      char **seq, int n, BioAlphabet &alphabet, 
      int mer_len, KmerGraph &kmer_freq, std::set<int> &uncorrectable,
      int max_iter = 100
  );
 private:
};

#endif
