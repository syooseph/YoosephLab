#ifndef _LOADER_H_
#define _LOADER_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cstring>
#include <cassert>
#include <list>
#include <unordered_map>

#include "bio_alphabet.h"

#ifndef MAXSEQLEN
#define MAXSEQLEN 1000000
#endif

class Loader  {
 public:
  Loader();
  ~Loader();
  // return the number of sequences in the FASTA file
  // file_name: the name of the file to be read
  int CountFastaNumSeqs(const char *file_name);
  // loading the FASTA file 
  // return the number of sequences
  // alphabet: the alphabet for the file
  // file_name: the name of the file to be read
  // header: the two-dimensional array that stores the header information
  // seq: the two-dimensional array that stores the sequence information
  int LoadFasta(BioAlphabet &alphabet, const char *file_name, char **header, char **seq);
 private:
  // checking for special characters in the string
  // alphabet: the alphabet for the string
  // sseq: the string to be checked
  // freq_cutoff: the portion of non-standard char cutoff to drop the string (return FALSE)
  //    otherwise modify the string and returns TRUE
  bool CheckSpecialChar(BioAlphabet &alphabet, std::string &sseq, float freq_cutoff = 0.9);
  void RecordSequence(
    char **header, char **seq, 
    std::string &single_header, std::string &single_seq, 
    int index
  );
};

#endif
