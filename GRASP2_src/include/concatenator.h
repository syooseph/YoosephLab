#ifndef _CONCATENATOR_H_
#define _CONCATENATOR_H_

#include <iostream>
#include <string>
#include <cstring>
#include <cassert>
#include <list>

typedef int64_t POSIDX;
typedef int32_t POSINT;

// the deliminator must be lexicographically the smallest
static char DELIM = (char) 1;

class Concatenator  {
 public:
  explicit Concatenator(char** const seq, const int n, std::string &concat_seq);
  explicit Concatenator(std::list<std::string> &seq, std::string &concat_seq);
  ~Concatenator() {}
};

#endif
