#include <iostream>

#include "../include/kmer.h"

using namespace std;

Kmer::Kmer(KmerUnitcoder &unitcoder, const std::string &seq)  {
  if(seq.length() <= 0)  {
    kmer_array_ = NULL;
    initialized_ = false;
    len_ = size_ = 0;
    return;
  }
  // record signature and the length
  this->signature_ = unitcoder.signature_;
  this->len_ = seq.length();
  // compute number of unit k-mer required for this k-mer
  this->size_ = (int) (len_ / unitcoder.mer_len_) + 1;
  int trailing_len = unitcoder.mer_len_ - len_ % unitcoder.mer_len_;
  // allocate space
  kmer_array_ = new KmerUnitType[size_];
  initialized_ = true;
  // encode
  int i;
  for(i = 0; i < size_ - 1; ++ i) {
    kmer_array_[i] = unitcoder.Encode(&(seq.c_str()[i * unitcoder.mer_len_]));
  }
  kmer_array_[i] = unitcoder.Encode(&(seq.c_str()[i * unitcoder.mer_len_ - trailing_len]));
  return;
}

Kmer::~Kmer() {
  if(initialized_) delete [] kmer_array_;
}

Kmer& Kmer::operator= (const Kmer &km)  {
  if(km.initialized_)  {
    if(this->initialized_ && this->size_ != km.size_)  {
      // different kmer should not be assigned
      cout << "Kmer::Error: assignment with k-mers with different lengths!" << endl;
      exit(1);
    } else if(!this->initialized_) {
      // otherwise if the destination is not initialized
      this->kmer_array_ = new KmerUnitType[km.size_];
    }
    for(int i = 0; i < km.size_; ++ i) {
      this->kmer_array_[i] = km.kmer_array_[i];
    }
  }
  this->signature_ = km.signature_;
  this->len_ = km.len_;
  this->size_ = km.size_;
  return *this;
}

std::string Kmer::Decode(KmerUnitcoder &unitcoder)  {
  if(this->signature_ != unitcoder.signature_)  {
    cout << "Error: attempting to decode k-mer failed, the k-mer is encoded and decoded using different alphabet and length settings." << endl;
    exit(1);
  }
  string mer = "";
  int trailing_len = unitcoder.mer_len_ - len_ % unitcoder.mer_len_;
  int i;
  for(i = 0; i < size_ - 1; ++ i) {
    mer += unitcoder.Decode(kmer_array_[i]);
  }
  mer += unitcoder.Decode(kmer_array_[i]).substr(trailing_len);
  return mer;
}

char Kmer::HeadChar(KmerUnitcoder &unitcoder) {
  if(this->signature_ != unitcoder.signature_)  {
    cout << "Error: attempting to decode k-mer failed, the k-mer is encoded and decoded using different alphabet and length settings." << endl;
    exit(1);
  }
  return unitcoder.HeadChar(kmer_array_[0]);
}

char Kmer::TailChar(KmerUnitcoder &unitcoder) {
  if(this->signature_ != unitcoder.signature_)  {
    cout << "Error: attempting to decode k-mer failed, the k-mer is encoded and decoded using different alphabet and length settings." << endl;
    exit(1);
  }
  return unitcoder.TailChar(kmer_array_[(int) len_ / unitcoder.mer_len_]);
}
