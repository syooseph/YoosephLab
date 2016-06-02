#ifndef _KMER_H_
#define _KMER_H_

#include <cstdlib>

#include "kmer_unitcoder.h"

class Kmer {
 public:
  explicit Kmer() :
    kmer_array_(NULL),
    initialized_(false),
    len_(0),
    size_(0)
  {}
  explicit Kmer(KmerUnitcoder &unitcoder, const std::string &seq);  
  ~Kmer();
  std::string Decode(KmerUnitcoder &unitcoder);
  // Get the first character of the k-mer
  char HeadChar(KmerUnitcoder &unitcoder);
  // Get the last character of the k-mer
  char TailChar(KmerUnitcoder &unitcoder);
  // overloading operator "=="
  bool operator==(const Kmer &km)  {
    for(int i = 0; i < size_; ++ i)
      if(this->kmer_array_[i] != km.kmer_array_[i]) return false;
    return true;
  }
  Kmer& operator=(const Kmer &km);
  // friend class definition 
  friend class DeBruijnGraph;
 protected:
  bool initialized_;
  KmerUnitType *kmer_array_;  
  KmerUnitType signature_;
  // the length of the k-mer
  uint8_t len_;
  // the size of the kmer_array_
  uint8_t size_;
};

#endif
