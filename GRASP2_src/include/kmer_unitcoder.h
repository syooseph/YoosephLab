#ifndef _KMER_UNITCODER_H_
#define _KMER_UNITCODER_H_

#include <iostream>
#include <map>
#include <string>
#include <cmath>
#include <cstdlib>
#include <bitset>
#include "bio_alphabet.h"

// using 32-bit unitcoder
typedef uint32_t KmerUnitType;

class KmerUnitcoder{
 public:
  // Object initializer
  // a: the alphabet of the sequence
  // l: the length of the k-mer
  explicit KmerUnitcoder(const BioAlphabet &a, const int l);
  ~KmerUnitcoder();
  // Encode a short string
  // s: the address of the begining of the string
  KmerUnitType Encode(const char *s);
  int EncodeInt(const char *s);
  int EncodeIntRight(const int num, const char c);
  // Decode an encoded string
  // ek: encoded k-mer
  std::string Decode(const KmerUnitType ek);
  // Encode the new k-mer constructed by concaternating the suffix (k-1)-mer and c
  KmerUnitType RightExt(const KmerUnitType ek, const char c);
  // Encode the new k-mer constructed by concaternating c and the prefix (k-1)-mer
  KmerUnitType LeftExt(const KmerUnitType ek, const char c);
  // Get the head character in the unit k-mer
  char HeadChar(const KmerUnitType ek);
  // Get the tail character in the unit k-mer
  char TailChar(const KmerUnitType ek);
  // Print out the unit k-mer as bits
  void PrintBits(const KmerUnitType ek);
  // Get kmer-length
  inline int GetMerLen(void) {return mer_len_;}
  // Friend class definition
  friend class Kmer;
  friend class DeBruijnGraph;
 protected:
  int mer_len_;
  int num_bit_;
  BioAlphabet alphabet_;
  KmerUnitType mask_;
  KmerUnitType signature_;    // the signature to make sure consistent unitcoder call
  void CkValid();
};

#endif
