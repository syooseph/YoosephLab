#include "../include/bwt.h"

using namespace std;


void BWT::BuildIndex(const BWTCHAR* text)  {
  bwt_ = new BWTCHAR [size_ + 1];
  int fm_size = (size_ / BWT_FM_GAP) + 1;
  int i, j;
  // allocate space for FM-index
  for(int i = 0; i < alphabet_.GetSize(); ++ i)  {
    fm_index_[alphabet_.GetInvCharMap(i)] = new BWTIDX [fm_size];
    for(j = 0; j < fm_size; ++ j) fm_index_[alphabet_.GetInvCharMap(i)][j] = 0;
  }
  fm_index_[DELIM] = new BWTIDX [fm_size];
  for(j = 0; j < fm_size; ++ j) fm_index_[DELIM][j] = 0;
  // count frequencies
  for(i = 0; i < size_; ++ i) {
    char c;
    if(position_[i] == 0) c = DELIM;
    else c = text[position_[i] - 1];
    // pad on the "bwt_" array
    bwt_[i] = c;  
    if(i % BWT_FM_GAP == 0)  {
      // copy the frequency to the FM index
      int ix = i / BWT_FM_GAP;
      for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
        if(fm_index_.find(it->first) == fm_index_.end())  {
          cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
          exit(1);
        }
        // add the frequency to the FM-index
        fm_index_[it->first][ix] = it->second;
      }
    }
    // add frequency
    ++ char_freq_[c];
  }
  bwt_[i] = '\0';
  // handle spacial case
  if(i % BWT_FM_GAP == 0) {
    // copy the frequency to the FM index
    int ix = i / BWT_FM_GAP;
    for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
      if(fm_index_.find(it->first) == fm_index_.end())  {
        cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
        exit(1);
      }
      // add the frequency to the FM-index
      fm_index_[it->first][ix] = it->second;
    }
  }  
  
  
  // compute accumulated index
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(i = 0; i < alphabet_.GetSize(); ++ i) {
    acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
  }

  // check if character that do not found in the alphabet exists
  for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
    if(it->first != DELIM && !alphabet_.IsValid((char) it->first))  {
      cerr << "Error: BWT::BuildIndex: Unrecognized character found in the text/sequence: " << it->first << "; abort." << endl;
      exit(1);
    }
  }
  return;
}

void BWT::BuildIndexNoPos(const BWTCHAR* text)  {
  int fm_size = (size_ / BWT_FM_GAP) + 1;
  int i, j;
  // allocate space for FM-index
  for(int i = 0; i < alphabet_.GetSize(); ++ i)  {
    fm_index_[alphabet_.GetInvCharMap(i)] = new BWTIDX [fm_size];
    for(j = 0; j < fm_size; ++ j) fm_index_[alphabet_.GetInvCharMap(i)][j] = 0;
  }
  fm_index_[DELIM] = new BWTIDX [fm_size];
  for(j = 0; j < fm_size; ++ j) fm_index_[DELIM][j] = 0;
  // count frequencies
  for(i = 0; i < size_; ++ i) { 
    if(i % BWT_FM_GAP == 0)  {
      // copy the frequency to the FM index
      int ix = i / BWT_FM_GAP;
      for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
        if(fm_index_.find(it->first) == fm_index_.end())  {
          cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
          exit(1);
        }
        // add the frequency to the FM-index
        fm_index_[it->first][ix] = it->second;
      }
    }
    // add frequency
    ++ char_freq_[bwt_[i]];
  }
  
  // handle spacial case
  if(i % BWT_FM_GAP == 0) {
    // copy the frequency to the FM index
    int ix = i / BWT_FM_GAP;
    for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
      if(fm_index_.find(it->first) == fm_index_.end())  {
        cerr << "Error: BWT::BuildIndex: Unrecognized character in text/sequence: " << it->first << "; abort." << endl;
        exit(1);
      }
      // add the frequency to the FM-index
      fm_index_[it->first][ix] = it->second;
    }
  }  
  
  
  // compute accumulated index
  acc_freq_ = new BWTIDX [alphabet_.GetSize() + 2]; // note the trailing character '$'
  // the begin of '$' in the array is 0
  acc_freq_[0] = 0; acc_freq_[1] = char_freq_[DELIM];
  for(i = 0; i < alphabet_.GetSize(); ++ i) {
    acc_freq_[i + 2] = acc_freq_[i + 1] + char_freq_[alphabet_.GetInvCharMap(i)];
  }

  // check if character that do not found in the alphabet exists
  for(auto it = char_freq_.begin(); it != char_freq_.end(); ++ it) {
    if(it->first != DELIM && !alphabet_.IsValid((char) it->first))  {
      cerr << "Error: BWT::BuildIndex: Unrecognized character found in the text/sequence: " << it->first << "; abort." << endl;
      exit(1);
    }
  }
  return;
}

void BWT::Construct(const BioAlphabet &alphabet, const char *text) {
  alphabet_ = alphabet;
  size_ = strlen(text);
  // build suffix array 
  position_ = new BWTIDX [size_]; 
  int divsufsort_tag = divsufsort64((BWTCHAR*) text, position_, size_);
#ifdef DEBUG
  cout << "Orders:  " << endl;
  for(int i = 0; i < size_; ++ i) {
    cout << position_[i] << " ";
  }
  cout << endl;
#endif
  // TODO: Handle sorting errors
  assert(divsufsort_tag == 0);
  // build indexes for BWT
  BuildIndex((BWTCHAR*) text);   
  build_success_ = true;
  return;
}

void BWT::ConstructNoPos(const BioAlphabet &alphabet, const char *text) {
  alphabet_ = alphabet;
  size_ = strlen(text);
  // build suffix array 
  bwt_ = new BWTCHAR [size_ + 1];
  int divsufsort_tag = divbwt64((BWTCHAR*) text, bwt_, NULL, size_);
  bwt_[size_] = '\0';
  // TODO: Handle sorting errors
#ifdef DEBUG
  cout << bwt_ << endl;
#endif
  assert(divsufsort_tag >= 0);
  // build indexes for BWT
  BuildIndexNoPos((BWTCHAR*) text);
  build_success_ = true;
  return;
}

std::pair<BWTIDX, BWTIDX> BWT::UpdateRange(
  const char c, const std::pair<BWTIDX, BWTIDX> &range
)  {
  int c_id = alphabet_.GetCharMap(c);  // the ID of the character in the alphabet
  if(c != DELIM && c_id < 0)  {
    // make sure that the character is in the alphabet
    cerr << "Unknown character being searched against the database: " << c << ", abort." << endl;
    exit(1);
  }
  // calculate the first occurrence
  BWTIDX occbegin = CountOccurrence(c, range.first);
  BWTIDX occend = CountOccurrence(c, range.second);
  // compute the new range
  pair<BWTIDX, BWTIDX> update_range;
  update_range.first = acc_freq_[c_id + 1] + occbegin;
  update_range.second = acc_freq_[c_id + 1] + occend;
  return update_range;  
}

BWTIDX BWT::CountLexicoLess(const char c, const BWTIDX pos) {
  BWTIDX occ = CountOccurrence(DELIM, pos);
  //cout << "CountLexicoLess: " << occ << endl;
  for(int i = 0; i < alphabet_.GetCharMap(c); ++ i) { 
    occ += CountOccurrence(alphabet_.GetInvCharMap(i), pos);
    //cout << "CountLexicoLess: " << alphabet_.GetInvCharMap(i) << "  " << occ << endl;
  }
  //cout << "End of function CountLexicoLess" << endl;
  return occ;
}
