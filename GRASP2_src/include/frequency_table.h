#ifndef _FREQUENCY_TABLE_H_
#define _FREQUENCY_TABLE_H_

#include <iostream>
#include <unordered_map>

template <class K_Type>
class FrequencyTable  {
 public:
  FrequencyTable(){};
  ~FrequencyTable(){};
  
  void Increase(const K_Type &k);
  void Decrease(const K_Type &k);
  void SetFreq(const K_Type &k, int f);
  int GetFreq(const K_Type &k) const;
  int GetSize(void) const;
  void PrintAll(void);
  
 private:
  std::unordered_map<K_Type, int> frequency_;
};

template <class K_Type>
inline void FrequencyTable<K_Type>::Increase(const K_Type &k) {
  auto it = frequency_.find(k);
  if(it != frequency_.end()) ++ it->second;
  else frequency_[k] = 1;
}

template <class K_Type>
inline void FrequencyTable<K_Type>::Decrease(const K_Type &k) {
  auto it = frequency_.find(k);
  if(it != frequency_.end()) {
    -- it->second;
    if(it->second <= 0) frequency_.erase(it);
  }
}

template <class K_Type>
inline void FrequencyTable<K_Type>::SetFreq(const K_Type &k, int f) {
  auto it = frequency_.find(k);
  if(it != frequency_.end()) it->second = f;
  else frequency_[k] = f;
}

template <class K_Type>
inline int FrequencyTable<K_Type>::GetFreq(const K_Type &k) const {
  auto it = frequency_.find(k);
  if(it != frequency_.end()) return it->second;
  else return 0;
}

template <class K_Type>
inline int FrequencyTable<K_Type>::GetSize(void) const {
  return frequency_.size();
}

template <class K_Type>
void FrequencyTable<K_Type>::PrintAll(void) {
  for(auto it = frequency_.begin(); it != frequency_.end(); ++ it)
    std::cout << it->first << "  " << it->second << std::endl;
  return;
}

#endif
