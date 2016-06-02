#ifndef _KMER_FILTERING_H_
#define _KMER_FILTERING_H_

#include "bio_alphabet.h"
#include "kmer_unitcoder.h"

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <unordered_map>
#include <stack>
#include <bitset>

const static int num_filters = 8;
const static int filter_size = 256;
const static int primes[8] = {251, 241, 239, 233, 229, 227, 223, 211}; 
//const static int filter_size = 1024;
//const static int primes[8] = {1021, 1019, 1013, 1009, 997, 991, 983, 977}; 
//const static int filter_size = 2048;
//const static int primes[8] = {2039, 2029, 2027, 2017, 2011, 2003, 1999, 1977}; 

class KmerFiltering {
 public:
  explicit KmerFiltering(int mer_len)  {
    mer_len_ = mer_len; init_ = false;      
  }
  explicit KmerFiltering()  {init_ = false;}
  ~KmerFiltering() {}
  // Computes the bloom-filter for all kmers (with length "mer_len") in string "seq"
  // "filter_size" is the size (number of bits) of the bloom filter
  void KmerBloomFilter(BioAlphabet &alphabet, const std::string &seq);
  // Computes the bloom-filter for all kmers (with length "mer_len") in the set "kmers"
  // "filter_size" is the size (number of bits) of the bloom filter
  void KmerBloomFilter(BioAlphabet &alphabet, std::vector<std::string> &kmers);
  // Compute which kmers are clustered together in the hash table
  // the first dimension of the km
  void GetKmerHashCluster(
      BioAlphabet &alphabet, std::vector<std::string> &kmers, 
      std::vector<std::vector<std::string> > &kmer_cluster    
  );
  // Convert a filter (an bool array) to a string
  std::string FilterToString(bool *p_filter);
  // Decode a coded string back to a filter
  void StringToFilter(std::string &str_coded, bool *p_filter);
  // check if all filters are proper loaded/prepared
  bool CheckFilterReady(void);
  // check if the encoded kmer unit exists in the set of kmers where the filters were built
  inline bool IsKmerUnitPresent(const KmerUnitType &k_unit)  {
    int id = 0;
    for(int i = 0; i < num_filters; ++ i) {
      if(!filters_[id + k_unit % primes[i]]) return false;
      id += filter_size;    
    }
    return true;
  }
  inline bool IsKmerUnitPresent(std::vector<int> &pos)  {
    for(int i = 0; i < pos.size(); ++ i) {
      if(!filters_[pos[i]]) return false;
    }
    return true;
  }
  // perform bit-wise comparison with another bloom filter and perform bit-wise AND
  // operation, which finds the positions that were set to 1 in both filters. The positions
  // are stored in pos
  void FindCoExistPos(KmerFiltering &bmfilter, std::vector<int> &pos);
  // find out which entries (as ID) in the set has a kmer being hashed at each position
  // the first dimension of "entries" is all sites in the filter (i.e. num_filters * filter_size)
  // the second dimension is a set of IDs in "bmfilters"
  void FindCoExistPosSet(
      std::vector<KmerFiltering> &bmfilters, std::vector<std::vector<int> > &entries
  );
  void ComputeSingleKmerHash(BioAlphabet &alphabet, const std::string &kmer, std::vector<int> &kvalues);
  // copy the filters into "filters"
  void GenFilterString(std::vector<std::string> &filter_str);
  int GetNumFilters(void) {return num_filters;}
  int GetFilterSize(void) {return filter_size;}
  void LoadFromStringSet(std::vector<std::string> &coded_filter, const int begin, const int end);
  inline void Purge(void)  {if(init_) delete [] filters_; init_ = false;}
  void PrintRawFilter(void) {
    for(int i = 0; i < filter_size * num_filters; ++ i) std::cout << filters_[i];
    std::cout << std::endl;
    return;
  }
  void PrintTureBitLocations()  {
    for(int i = 0; i < filter_size * num_filters; ++ i) {
      if(filters_[i]) std::cout << ">>>:  " << i << std::endl;
    }
    return;
  }
 private:
  int mer_len_; bool init_;
  bool *filters_;
};

#endif
