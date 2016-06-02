#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <string>
#include <list>

#ifndef _REDUCED_ALPHABET_
#define _REDUCED_ALPHABET_

// a list of supported alphabets
enum Alphabet {ALL20, DSSP5, DSSP10, GBMR4, GBMR10, HSDM5, SDM6, MURPHY5, MURPHY10, TD5, TD10};

class ReducedAlphabet {
 public:
  ReducedAlphabet(enum Alphabet a);
  ~ReducedAlphabet(void);
  // convert from original alphabet to the reduced alphabet
  char Convert(char c, bool& is_standard); 
  std::string Convert(std::string s, bool& is_standard); 
 protected:
  // data begin
  enum Alphabet alphabet_in_use_;
  std::map<char, int> alphabet_map_;
  // methods begin
  void InitALL20(void);
  void InitDSSP5(void);
  void InitDSSP10(void);
  void InitGBMR4(void);
  void InitGBMR10(void);
  void InitHSDM5(void);
  void InitSDM6(void);
  void InitMURPHY5(void);
  void InitMURPHY10(void);
  void InitTD5(void);
  void InitTD10(void); 
};

#endif
