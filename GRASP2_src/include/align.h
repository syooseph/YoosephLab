#ifndef _ALIGN_H_
#define _ALIGN_H_

#include <iostream>
#include <string>
#include <vector>

// define the maximum edit distance
static const int dist_max = 999999;

class Align {
 public:
  explicit Align() {};
  ~Align() {};
  // golbal alignment with edit distance
  int GlobalAlnEdd(const std::string &seqA, const std::string &seqB, int band_size = 3);
  // Sequence A must be global at both begin and end (GG)
  // Sequence B must be global at begin, but can be local and end (GL)
  int GGGLAlnEdd(const std::string &seqA, const std::string &seqB, int band_size = 3);
};


int Align::GlobalAlnEdd(const std::string &seqA, const std::string &seqB, int band_size) {
  if(seqA.length() == 0) return seqB.length();
  else if(seqB.length() == 0) return seqA.length();
  // initialize the table to all '0's
  int l_band = 2 * band_size + 1;
  std::vector<int> T(l_band * (seqA.length() + 1), dist_max);
  // initialize
  T[band_size] = 0;
  for(int i = 0; i < band_size; ++ i) {T[i] = dist_max; T[band_size + i + 1] = i + 1;}
  // fill the table
  int pivot = 0;
  int min_dist = dist_max;
  for(int i = 1; i <= seqA.length(); ++ i) { 
    // note i <= seqA.length() because the first row represents deletion
    if(i > seqB.length() + band_size)  break;
    pivot += l_band;
    for(int j = 0; j < l_band; ++ j) {
      //std::cout << pivot << " " << j << std::endl;
      int j_idx = i + j - band_size;  // the actual position of j in the hypothesized 2D table
      //std::cout << "Actual J: " << j_idx << std::endl;
      if(j_idx < 0 || j_idx > seqB.length()) continue;
      // evaluate a match
      if(j_idx > 0 && T[pivot + j] > T[pivot + j - l_band] + (seqA[i - 1] != seqB[j_idx - 1]))  {  
        //std::cout << "Diagnoal: " << T[pivot + j - l_band] << std::endl;
        T[pivot + j] = T[pivot + j - l_band] + (seqA[i - 1] != seqB[j_idx - 1]);
      }
      if(j > 0 && T[pivot + j] > T[pivot + j - 1] + 1)  {
        //std::cout << "Left: " << T[pivot + j - 1] << std::endl;
        T[pivot + j] = T[pivot + j - 1] + 1;
      }
      //std::cout << "Upper score: " << T[pivot + j] << "  " << T[pivot + j - l_band + 1] << std::endl;
      if(j < l_band - 1 && T[pivot + j] > T[pivot + j - l_band + 1] + 1)  {
        //std::cout << "Upper: " << T[pivot + j - l_band + 1] << std::endl;
        T[pivot + j] = T[pivot + j - l_band + 1] + 1;
      }
      int ext_dist = seqA.length() - i > seqB.length() - j_idx ? seqA.length() - i : seqB.length() - j_idx;
      if(min_dist > T[pivot + j] + ext_dist)
        min_dist = T[pivot + j] + ext_dist;
    }
  }
  
  // DEBUG
  //for(auto it = T.begin(); it != T.end(); ++ it)
  //  std::cout << *it << "\t";
  //std::cout << std::endl;
  return min_dist;
}


int Align::GGGLAlnEdd(const std::string &seqA, const std::string &seqB, int band_size)  {
  if(seqA.length() == 0) return seqB.length();
  else if(seqB.length() == 0) return seqA.length();
  // initialize the table to all '0's
  int l_band = 2 * band_size + 1;
  std::vector<int> T(l_band * (seqA.length() + 1), dist_max);
  // initialize
  T[band_size] = 0;
  for(int i = 0; i < band_size; ++ i) {T[i] = dist_max; T[band_size + i + 1] = i + 1;}
  // fill the table
  int pivot = 0;
  int min_dist = dist_max;
  for(int i = 1; i <= seqA.length(); ++ i) { 
    // note i <= seqA.length() because the first row represents deletion
    if(i > seqB.length() + band_size)  break;
    pivot += l_band;
    for(int j = 0; j < l_band; ++ j) {
      //std::cout << pivot << " " << j << std::endl;
      int j_idx = i + j - band_size;  // the actual position of j in the hypothesized 2D table
      //std::cout << "Actual J: " << j_idx << std::endl;
      if(j_idx < 0 || j_idx > seqB.length()) continue;
      // evaluate a match
      if(j_idx > 0 && T[pivot + j] > T[pivot + j - l_band] + (seqA[i - 1] != seqB[j_idx - 1]))  {  
        //std::cout << "Diagnoal: " << T[pivot + j - l_band] << std::endl;
        T[pivot + j] = T[pivot + j - l_band] + (seqA[i - 1] != seqB[j_idx - 1]);
      }
      if(j > 0 && T[pivot + j] > T[pivot + j - 1] + 1)  {
        //std::cout << "Left: " << T[pivot + j - 1] << std::endl;
        T[pivot + j] = T[pivot + j - 1] + 1;
      }
      //std::cout << "Upper score: " << T[pivot + j] << "  " << T[pivot + j - l_band + 1] << std::endl;
      if(j < l_band - 1 && T[pivot + j] > T[pivot + j - l_band + 1] + 1)  {
        //std::cout << "Upper: " << T[pivot + j - l_band + 1] << std::endl;
        T[pivot + j] = T[pivot + j - l_band + 1] + 1;
      }
      int ext_dist = seqA.length() - i; // we don't count the cost of seqB
      if(min_dist > T[pivot + j] + ext_dist)
        min_dist = T[pivot + j] + ext_dist;
    }
  }
  
  // DEBUG
  //for(auto it = T.begin(); it != T.end(); ++ it)
  //  std::cout << *it << "\t";
  //std::cout << std::endl;
  return min_dist;
}
#endif








