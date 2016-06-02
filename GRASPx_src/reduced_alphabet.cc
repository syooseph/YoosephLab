#include "reduced_alphabet.h"

using namespace std;

ReducedAlphabet::ReducedAlphabet(enum Alphabet a) {
  alphabet_in_use_ = a;
  switch(a)  {
    case ALL20: InitALL20(); break;
    case DSSP5: InitDSSP5(); break;
    case DSSP10: InitDSSP10(); break;
    case GBMR4: InitGBMR4(); break;
    case GBMR10: InitGBMR10(); break;
    case HSDM5: InitHSDM5(); break;
    case SDM6: InitSDM6(); break;
    case MURPHY5: InitMURPHY5(); break;
    case MURPHY10: InitMURPHY10(); break;
    case TD5: InitTD5(); break;
    case TD10: InitTD10(); break;
    default: 
      cout << "Error: ReducedAlphabet::ReducedAlphabet alphabet currently not supported" << endl; 
      exit(0);
  }
  return;
}

ReducedAlphabet::~ReducedAlphabet(void)  {
  return;
}

char ReducedAlphabet::Convert(char c, bool& is_standard)  {
  is_standard = true;
  if(alphabet_map_.find(c) != alphabet_map_.end())  {
    // 65 is added to convert the reduced alphabet to readable characters
    return char (65 + alphabet_map_[c]);
  }  else  {
    is_standard = false;
    return 'X';
  }
}

std::string ReducedAlphabet::Convert(std::string s, bool& is_standard)  {
  is_standard = true;
  string n = "";
  for(auto it = s.begin(); it != s.end(); ++ it) {
    bool is_std_loc;
    n += Convert(*it, is_std_loc);
    // mark as non-standard is any of the character is non-standard
    if(!is_std_loc)  {
      is_standard = false;
    }
  }
  return n;
}

void ReducedAlphabet::InitALL20(void)  {
  alphabet_map_ = {
    {'P',  0}, {'G',  1}, {'E',  2}, {'K',  3}, {'R',  4}, 
    {'Q',  5}, {'D',  6}, {'S',  7}, {'N',  8}, {'T',  9}, 
    {'H', 10}, {'C', 11}, {'I', 12}, {'V', 13}, {'W', 14}, 
    {'Y', 15}, {'F', 16}, {'A', 17}, {'L', 18}, {'M', 19}
  };
  return;
}

void ReducedAlphabet::InitDSSP5(void)  {
  alphabet_map_ = {
    {'P',  4}, {'G',  4}, {'E',  0}, {'K',  0}, {'R',  0}, 
    {'Q',  0}, {'D',  3}, {'S',  2}, {'N',  3}, {'T',  2}, 
    {'H',  0}, {'C',  2}, {'I',  1}, {'V',  1}, {'W',  1}, 
    {'Y',  1}, {'F',  1}, {'A',  0}, {'L',  1}, {'M',  1}
  };
  return;
}

void ReducedAlphabet::InitDSSP10(void)  {
  alphabet_map_ = {
    {'P',  9}, {'G',  9}, {'E',  0}, {'K',  0}, {'R',  0}, 
    {'Q',  0}, {'D',  8}, {'S',  8}, {'N',  8}, {'T',  6}, 
    {'H',  6}, {'C',  7}, {'I',  1}, {'V',  1}, {'W',  5}, 
    {'Y',  2}, {'F',  3}, {'A',  4}, {'L',  2}, {'M',  4}
  };
  return;
}

void ReducedAlphabet::InitGBMR4(void)  {
  alphabet_map_ = {
    {'P',  3}, {'G',  0}, {'E',  1}, {'K',  1}, {'R',  1}, 
    {'Q',  1}, {'D',  1}, {'S',  1}, {'N',  1}, {'T',  1}, 
    {'H',  2}, {'C',  2}, {'I',  2}, {'V',  2}, {'W',  2}, 
    {'Y',  2}, {'F',  2}, {'A',  1}, {'L',  2}, {'M',  2}
  };
  return;
}

void ReducedAlphabet::InitGBMR10(void)  {
  alphabet_map_ = {
    {'P',  9}, {'G',  0}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  1}, {'S',  8}, {'N',  2}, {'T',  7}, 
    {'H',  5}, {'C',  6}, {'I',  3}, {'V',  3}, {'W',  3}, 
    {'Y',  4}, {'F',  3}, {'A',  3}, {'L',  3}, {'M',  3}
  };
  return;
}

void ReducedAlphabet::InitHSDM5(void)  {
  alphabet_map_ = {
    {'P',  3}, {'G',  3}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  3}, 
    {'H',  4}, {'C',  2}, {'I',  0}, {'V',  0}, {'W',  1}, 
    {'Y',  0}, {'F',  0}, {'A',  3}, {'L',  0}, {'M',  0}
  };
  return;
}

void ReducedAlphabet::InitSDM6(void)  {
  alphabet_map_ = {
    {'P',  5}, {'G',  3}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  3}, 
    {'H',  4}, {'C',  1}, {'I',  0}, {'V',  0}, {'W',  2}, 
    {'Y',  0}, {'F',  0}, {'A',  3}, {'L',  0}, {'M',  0}
  };
  return;
}

void ReducedAlphabet::InitMURPHY5(void)  {
  alphabet_map_ = {
    {'P',  1}, {'G',  1}, {'E',  3}, {'K',  4}, {'R',  4}, 
    {'Q',  3}, {'D',  3}, {'S',  1}, {'N',  3}, {'T',  1}, 
    {'H',  4}, {'C',  0}, {'I',  0}, {'V',  0}, {'W',  2}, 
    {'Y',  2}, {'F',  2}, {'A',  1}, {'L',  0}, {'M',  0}
  };
  return;
}

void ReducedAlphabet::InitMURPHY10(void)  {
  alphabet_map_ = {
    {'P',  8}, {'G',  4}, {'E',  2}, {'K',  1}, {'R',  1}, 
    {'Q',  2}, {'D',  2}, {'S',  9}, {'N',  2}, {'T',  9}, 
    {'H',  5}, {'C',  3}, {'I',  6}, {'V',  6}, {'W',  7}, 
    {'Y',  7}, {'F',  7}, {'A',  0}, {'L',  6}, {'M',  6}
  };
  return;
}

void ReducedAlphabet::InitTD5(void)  {
  alphabet_map_ = {
    {'P',  0}, {'G',  0}, {'E',  1}, {'K',  1}, {'R',  1}, 
    {'Q',  1}, {'D',  2}, {'S',  2}, {'N',  2}, {'T',  2}, 
    {'H',  2}, {'C',  2}, {'I',  3}, {'V',  3}, {'W',  3}, 
    {'Y',  3}, {'F',  3}, {'A',  4}, {'L',  4}, {'M',  4}
  };
  return;
}

void ReducedAlphabet::InitTD10(void)  {
  alphabet_map_ = {
    {'P',  0}, {'G',  1}, {'E',  2}, {'K',  2}, {'R',  2}, 
    {'Q',  2}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  4}, 
    {'H',  5}, {'C',  5}, {'I',  6}, {'V',  6}, {'W',  7}, 
    {'Y',  7}, {'F',  7}, {'A',  8}, {'L',  9}, {'M',  9}
  };
  return;
}

