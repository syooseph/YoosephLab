#ifndef _BIO_ALPHABET_H_
#define _BIO_ALPHABET_H_

#include <iostream>
#include <map>
#include <cstdlib>
#include <cstring>
#include <tuple>
#include <time.h>
#include <vector>

enum BioSequence {PROT, DNA, RNA};

class BioAlphabet{
 public:
  explicit BioAlphabet() {}
  explicit BioAlphabet(const BioSequence s);
  ~BioAlphabet() {}
  // returns a random character within the alphabet
  char RandomChar();
  bool CheckSuccess();
  bool IsValid(const char c);
  inline int GetSize(void) {return alphabet_size_;}
  inline int GetCharMap(const char c) {
    return char_map_[(int) c];
  }
  inline char GetInvCharMap(const int n) {
    return (char) inv_char_map_[n];
  }
  
  BioAlphabet& operator=(const BioAlphabet &alpha_source) {
    this->seq_type_ = alpha_source.seq_type_;
    this->alphabet_size_ = alpha_source.alphabet_size_;
    this->char_map_ = alpha_source.char_map_;
    this->inv_char_map_ = alpha_source.inv_char_map_;
    return *this;
  }
  
  friend class KmerUnitcoder;
  friend class BWT;
  friend class BWTSearch;

 protected:
  void InitProt();
  void InitDNA();
  void InitRNA();
  std::vector<int> char_map_;
  std::vector<int> inv_char_map_;
  BioSequence seq_type_;
  int alphabet_size_;
  
};
#endif
