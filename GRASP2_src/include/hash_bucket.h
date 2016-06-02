#ifndef _HASH_BUCKET_H_
#define _HASH_BUCKET_H_

#include <unordered_map>
#include <list>

/*
  This class primarily handles the search, insert, and delete 
  of elements in a hash bucket. This class requires that the 
  V_Type class has the operated "==" overloaded.
*/

template<class K_Type, class V_Type>
class HashBucket  {
 public:
  explicit HashBucket(){}
  ~HashBucket(){}
  void Insert(const K_Type &key, const V_Type &value);
  void Delete(const K_Type &key, const V_Type &value);
  V_Type* Search(const K_Type &key, const V_Type &value);
 protected:
  std::unordered_map<K_Type, std::list<V_Type> > hash_;
};

template<class K_Type, class V_Type>
V_Type* HashBucket<K_Type, V_Type>::Search(
    const K_Type &key, const V_Type &value
) {
  auto it_unmap = hash_.find(key);
  if(it_unmap != hash_.end())  {
    // search each element in the bucket one-by-one
    auto it_v = it_unmap->second.begin();
    for(; it_v != it_unmap->second.end(); ++ it_v) {
      if(*it_v == value)  break;
    }
    if(it_v != it_unmap->second.end())  return &(*it_v);
  } 
  return NULL;
}

template<class K_Type, class V_Type>
inline void HashBucket<K_Type, V_Type>::Insert(const K_Type &key, const V_Type &value) {
  if(Search(key, value) == NULL) hash_[key].push_back(value);
}

template<class K_Type, class V_Type>
void HashBucket<K_Type, V_Type>::Delete(const K_Type &key, const V_Type &value) {
  auto it_unmap = hash_.find(key);
  int del_copies = 0;
  if(it_unmap != hash_.end())  {
    // search each element in the bucket one-by-one
    auto it_v = it_unmap->second.begin();    
    for(; it_v != it_unmap->second.end(); ++ it_v) {
      if(*it_v == value)  {++ del_copies; it_unmap->erase(it_v);}
    }
  }
  if(del_copies <= 0) 
      std::cout << "HashBucket::Warning: attemping to delete entry that does not exists!" << std::endl;
  return;
}

#endif
