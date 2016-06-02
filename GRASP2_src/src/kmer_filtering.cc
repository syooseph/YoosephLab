#include "../include/kmer_filtering.h"

using namespace std;

void KmerFiltering::ComputeSingleKmerHash(
    BioAlphabet &alphabet, const std::string &kmer, std::vector<int> &kvalues
)  {
  if(kmer.length() != mer_len_) {
    cout << "Error: KmerFiltering::ComputeSingleKmerHash: inconsistent kmer length, abort" << endl;
    exit(1);
  }  
  KmerUnitcoder kmer_coder(alphabet, mer_len_);
  KmerUnitType ck = kmer_coder.Encode(kmer.c_str());
  for(int k = 0; k < num_filters; ++ k) {
    kvalues.push_back(k * filter_size + ck % primes[k]);
  }
  return;
}

void KmerFiltering::KmerBloomFilter(BioAlphabet &alphabet, const std::string &seq) {
  if(seq.length() < mer_len_)  return;
  // construct k-mer unitcoder object
  filters_ = new bool [filter_size * num_filters];
  memset(filters_, false, filter_size * num_filters);
  init_ = true;
  int i, j, k, a = sizeof(KmerUnitType);
  KmerUnitcoder kmer_coder(alphabet, mer_len_);
  uint32_t tail_one = 1;
  KmerUnitType ck;
  for(i = 0; i <= seq.length() - mer_len_; ++ i) {
    if(i == 0)  ck = kmer_coder.Encode(seq.substr(0, mer_len_).c_str());
    else ck = kmer_coder.RightExt(ck, seq[i + mer_len_ - 1]);
    for(k = 0; k < num_filters; ++ k) {
      filters_[k * filter_size + ck % primes[k]] = true;
    }
  } 
  return;
}

void KmerFiltering::KmerBloomFilter(BioAlphabet &alphabet, std::vector<std::string> &kmers)  {
  if(kmers.size() <= 0)  return;
  // construct k-mer unitcoder object
  filters_ = new bool [filter_size * num_filters];
  memset(filters_, false, filter_size * num_filters);
  init_ = true;
  int i, j, k, a = sizeof(KmerUnitType);
  KmerUnitcoder kmer_coder(alphabet, mer_len_);
  uint32_t tail_one = 1;
  KmerUnitType ck;
  for(i = 0; i < kmers.size(); ++ i) {
    if(kmers[i].length() != mer_len_)  continue;
    ck = kmer_coder.Encode(kmers[i].c_str());
    for(k = 0; k < num_filters; ++ k) {
      filters_[k * filter_size + ck % primes[k]] = true;
    }
  }
  return;
}

void KmerFiltering::GetKmerHashCluster(
    BioAlphabet &alphabet, std::vector<std::string> &kmers, 
    std::vector<std::vector<std::string> > &kmer_cluster    
) {
  if(kmers.size() <= 0)  return;
  // construct k-mer unitcoder object
  kmer_cluster.resize(filter_size * num_filters);
  int i, j, k, a = sizeof(KmerUnitType);
  KmerUnitcoder kmer_coder(alphabet, mer_len_);
  uint32_t tail_one = 1;
  KmerUnitType ck;
  for(i = 0; i < kmers.size(); ++ i) {
    if(kmers[i].length() != mer_len_)  continue;
    ck = kmer_coder.Encode(kmers[i].c_str());
    for(k = 0; k < num_filters; ++ k) {
      kmer_cluster[k * filter_size + ck % primes[k]].push_back(kmers[i]);
    }
  }
  return;
}

void KmerFiltering::FindCoExistPos(KmerFiltering &bmfilter, std::vector<int> &pos)  {
  if(num_filters != bmfilter.GetNumFilters() || filter_size != bmfilter.GetFilterSize())  {
    cout << "Warning: KmerFiltering::FindCoExistPos: The two bloom filters have different parameters!!!" << endl;
    return;
  }
  // compare positions one by one;
  for(int i = 0; i < num_filters * filter_size; ++ i) {
    if(filters_[i] && bmfilter.filters_[i]) pos.push_back(i);
  }
  return;
}

void KmerFiltering::FindCoExistPosSet(
    std::vector<KmerFiltering> &bmfilters, std::vector<std::vector<int> > &entries
) {
  int i, j;
  int n = num_filters * filter_size;
  entries.resize(n);
  for(i = 0; i < bmfilters.size(); ++ i) {
    for(j = 0; j < n; ++ j) {
      if(bmfilters[i].filters_[j]) entries[j].push_back(i);
    }
  }
  return;
}

std::string KmerFiltering::FilterToString(bool *p_filter)  {
  int i, j;
  string out_s = "";
  for(i = 0; i * 8 < filter_size; ++ i) {
    uint8_t c = 0;
    for(j = 0; j < 8; ++ j) {
      c = c << 1;
      int id = i * 8 + j;
      uint8_t d = p_filter[id];
      c = c | d;
    }
    out_s += (char) c;
  }
  return out_s;
}

void KmerFiltering::GenFilterString(std::vector<std::string> &filter_str)  {
  for(int i = 0; i < num_filters; ++ i) {
    filter_str.push_back(FilterToString(&filters_[i * filter_size]));
  }
  return;
}

void KmerFiltering::StringToFilter(std::string &str_coded, bool *p_filter)  {
  // check string length
  int exs = filter_size / 8;
  if(str_coded.length() != exs) {
    cout << "Warning: KmerFiltering::StringToFilter: The length of the string is incompatible with the filter size!!!" << endl;  
    return;
  }
  // decode the string 
  uint8_t tail_one = 1;
  for(int i = 0; i < str_coded.length(); ++ i) {
    uint8_t c = (uint8_t) str_coded[i];
    int n = (i + 1) * 8 - 1;
    for(int j = 0; j < 8; ++ j) {
      uint8_t d = c & tail_one;
      p_filter[n - j] = (bool) d;
      c = c >> 1;
    }
  }
  return;
}

// check if all filters are proper loaded/prepared
bool KmerFiltering::CheckFilterReady(void)  {
  return init_;
}

void KmerFiltering::LoadFromStringSet(
    std::vector<std::string> &coded_filter, const int begin, const int end
) {
  if(end - begin + 1 != num_filters)  {
    cout << "Error: KmerFiltering::LoadFromStringSet: Inconsistent number of coded filter strings, abort." << endl;
    exit(1);
  }
  filters_ = new bool [filter_size * num_filters];
  memset(filters_, false, filter_size * num_filters);
  init_ = true;
  int i, j;
  for(i = begin, j = 0; i <= end; ++ i, ++ j) {
    StringToFilter(coded_filter[i], &filters_[j * filter_size]);
  }
  return;
}
