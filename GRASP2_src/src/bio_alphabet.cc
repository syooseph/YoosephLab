#include "../include/bio_alphabet.h"

using namespace std;

BioAlphabet::BioAlphabet(const BioSequence s) {
  seq_type_ = s;
  char_map_.assign(256, -1);
  inv_char_map_.assign(256, -1);
  srand(time(NULL));
  switch (s)  {
    case PROT: InitProt();  break;
    case DNA: InitDNA();  break;
    case RNA: InitRNA();  break;
    default: std::cout << "Unknown BIOSEQUENCE: the program only supports PROT, DNA, or RNA." << std::endl; break;
  }
  return;
}


void BioAlphabet::InitProt(void) {
  alphabet_size_ = 25;
  char_map_[(int) '*'] = 0; inv_char_map_[0] = (int) '*';
  char_map_[(int) 'A'] = 1; inv_char_map_[1] = (int) 'A';
  char_map_[(int) 'B'] = 2; inv_char_map_[2] = (int) 'B';
  char_map_[(int) 'C'] = 3; inv_char_map_[3] = (int) 'C';
  char_map_[(int) 'D'] = 4; inv_char_map_[4] = (int) 'D';
  char_map_[(int) 'E'] = 5; inv_char_map_[5] = (int) 'E';
  char_map_[(int) 'F'] = 6; inv_char_map_[6] = (int) 'F';  
  char_map_[(int) 'G'] = 7; inv_char_map_[7] = (int) 'G';
  char_map_[(int) 'H'] = 8; inv_char_map_[8] = (int) 'H';
  char_map_[(int) 'I'] = 9; inv_char_map_[9] = (int) 'I';
  char_map_[(int) 'J'] = 10; inv_char_map_[10] = (int) 'J';
  char_map_[(int) 'K'] = 11; inv_char_map_[11] = (int) 'K';  
  char_map_[(int) 'L'] = 12; inv_char_map_[12] = (int) 'L';
  char_map_[(int) 'M'] = 13; inv_char_map_[13] = (int) 'M';
  char_map_[(int) 'N'] = 14; inv_char_map_[14] = (int) 'N';
  char_map_[(int) 'P'] = 15; inv_char_map_[15] = (int) 'P';
  char_map_[(int) 'Q'] = 16; inv_char_map_[16] = (int) 'Q';
  char_map_[(int) 'R'] = 17; inv_char_map_[17] = (int) 'R';
  char_map_[(int) 'S'] = 18; inv_char_map_[18] = (int) 'S';
  char_map_[(int) 'T'] = 19; inv_char_map_[19] = (int) 'T';
  char_map_[(int) 'V'] = 20; inv_char_map_[20] = (int) 'V';
  char_map_[(int) 'W'] = 21; inv_char_map_[21] = (int) 'W';
  char_map_[(int) 'X'] = 22; inv_char_map_[22] = (int) 'X';
  char_map_[(int) 'Y'] = 23; inv_char_map_[23] = (int) 'Y';
  char_map_[(int) 'Z'] = 24; inv_char_map_[24] = (int) 'Z';
  
  return;
}

void BioAlphabet::InitDNA(void) {
  alphabet_size_ = 4;
  char_map_[(int) 'A'] = 0; inv_char_map_[0] = (int) 'A';
  char_map_[(int) 'C'] = 1; inv_char_map_[1] = (int) 'C';
  char_map_[(int) 'G'] = 2; inv_char_map_[2] = (int) 'G';
  char_map_[(int) 'T'] = 3; inv_char_map_[3] = (int) 'T';
  return;
}

void BioAlphabet::InitRNA(void) {
  alphabet_size_ = 4;
  char_map_[(int) 'A'] = 0; inv_char_map_[0] = (int) 'A';
  char_map_[(int) 'C'] = 1; inv_char_map_[1] = (int) 'C';
  char_map_[(int) 'G'] = 2; inv_char_map_[2] = (int) 'G';
  char_map_[(int) 'U'] = 3; inv_char_map_[3] = (int) 'U';
  return;
}

bool BioAlphabet::IsValid(const char c) {
  return (char_map_[(int) c] != -1);
}

char BioAlphabet::RandomChar()  {
  int pos = rand() % alphabet_size_;
  return (char) inv_char_map_[pos];
}
