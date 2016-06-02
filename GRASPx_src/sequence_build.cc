#include "sequence_build.h"

using namespace std;

SequenceBuild::SequenceBuild(void)  {
  is_header_loaded_ = is_sequence_loaded_ = false;
  is_sfa_built_ = is_k_array_built_ = false; 
  is_size_counted_ = false;
  return;
}

SequenceBuild::SequenceBuild(std::string& file_name)  {
  is_header_loaded_ = is_sequence_loaded_ = is_sfa_built_ = false;
  is_size_counted_ = false;
  LoadSequences(file_name);
  return;
}

SequenceBuild::SequenceBuild(SequenceBuild &seq_obj)  {
  this->is_sfa_built_ = false;
  this->is_size_counted_ = false;
  this->num_seqs_ = seq_obj.num_seqs_;
  this->header_ = new char* [num_seqs_];
  this->sequence_ = new char* [num_seqs_];
  for(int i = 0; i < num_seqs_; ++ i) {
    if(seq_obj.is_header_loaded_)  {
      // copy the header
      int n = strlen(seq_obj.header_[i]) + 1;
      this->header_[i] = new char [n];
      memcpy(this->header_[i], seq_obj.header_[i], n);
    }
    if(seq_obj.is_sequence_loaded_)  {
      // copy the sequence
      int n = strlen(seq_obj.sequence_[i]) + 1;
      this->sequence_[i] = new char [n];
      memcpy(this->sequence_[i], seq_obj.sequence_[i], n);
    }
  }
  
  this->is_header_loaded_ = seq_obj.is_header_loaded_;
  this->is_sequence_loaded_ = seq_obj.is_sequence_loaded_;
  return;
}


SequenceBuild::~SequenceBuild(void) {
  if(is_sequence_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] sequence_[i];
    }
    delete [] sequence_;
  }
  if(is_header_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] header_[i];
    }
    delete [] header_;
  }
  if(is_sfa_built_)  {
    delete suffix_array_;
  }
  return;
}

void SequenceBuild::DestructSFA(void)  {
  // clear suffix array and key array to save memory
  suffix_array_->clear();
  suffix_array_->purgeSA();
  suffix_array_->purgeDoc();
  suffix_array_->purgeLCP();
  suffix_array_->purgeMLCP();
  key_array_.clear();
  is_sfa_built_ = false;
  is_k_array_built_ = false;
}

void SequenceBuild::DestructSequences(void) {
  if(is_sequence_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] sequence_[i];
    }
    is_sequence_loaded_ = false;
  }
  if(is_header_loaded_)  {
    for(int i = 0; i < num_seqs_; ++ i) {
      delete [] header_[i];
    }
    is_header_loaded_ = false;
  }
  num_seqs_ = 0;
  return;
}

void SequenceBuild::PrintAllSeqs(void)  {
  cout << "Sequences in object:" << endl;
  for(int i = 0; i < num_seqs_; ++ i) {
    cout << sequence_[i] << endl;
  }
  return;
}

std::string SequenceBuild::GetFileStem(const std::string& path)  {
  // going backward untill the '\/' character
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    if(path[i] == '/')  {
      break;
    }
  }
  return path.substr(i + 1, path.length() - i - 1);
}

void SequenceBuild::DumpSFA(std::string& dir, std::string& file_stem) {
  string lcp_file = dir + "/" + file_stem + ".lcp";
  string mcp_file = dir + "/" + file_stem + ".mcp";
  string gsa_file = dir + "/" + file_stem + ".gsa";
  // ignoring SFA and DOC files
  suffix_array_->dump(
      "", "", lcp_file.c_str(), mcp_file.c_str(), gsa_file.c_str()
  );
  return;
}

void SequenceBuild::LoadSFA(std::string& dir, std::string& file_stem) {
  string lcp_file = dir + "/" + file_stem + ".lcp";
  string mcp_file = dir + "/" + file_stem + ".mcp";
  string gsa_file = dir + "/" + file_stem + ".gsa";
  // ignoring SFA and DOC files
  if(!boost::filesystem::exists(lcp_file) ||
     !boost::filesystem::exists(mcp_file) ||
     !boost::filesystem::exists(gsa_file))  {
    cerr << "Error: SequenceBuild::LoadSFA: Indexing files not found: have you run \'grasp-build\' first?" << endl;
    exit(0);
  } else  {
    // load the existing suffix array
    suffix_array_ = new GSA();
    suffix_array_->load(
        lcp_file.c_str(), mcp_file.c_str(), gsa_file.c_str()
    );
    suffix_array_->setSequences(sequence_);
    suffix_array_->setReadCount(num_seqs_);
    is_sfa_built_ = true;
  }
  return;
}

void SequenceBuild::CountDBSize(void) {
  if(is_size_counted_)  {
    return;
  }
  db_size_MB_ = 0;
  for(int i = 0; i < num_seqs_; ++ i) {
    db_size_MB_ += (double) strlen(sequence_[i]);
  }
  db_size_MB_ /= 1000000;
  is_size_counted_ = true;
  return;
}

double SequenceBuild::GetDBSizeInMegaBase(void)  {
  if(!is_size_counted_)  {
    this->CountDBSize();
  }
  return db_size_MB_;
}


// the sequences will be stored in global variables "header" and "sequence_"
void SequenceBuild::LoadSequences(std::string& file_name)  {
  vector<string> files_in;
  files_in.push_back(file_name);
  num_seqs_ = (unsigned int) seq::totalSequenceCount(files_in);
  header_ = new char* [num_seqs_];
  sequence_ = new char* [num_seqs_];
  seq::loadSequences(files_in, header_, sequence_, TAGSEQ);
  for(int i = 0; i < num_seqs_; ++ i) {
    int l = strlen(sequence_[i]);
    // chomp tailing non-characters
    while(!isalpha(sequence_[i][l - 1]))  {
      sequence_[i][l - 1] = '\0';
      -- l;
    }
  }
  is_header_loaded_ = is_sequence_loaded_ = true;
  return;
}

// copy sequence from one array to another
void SequenceBuild::CopySequences(int num_sequences, char **source, char **target) {
  target = new char* [num_sequences];
  for(int i = 0; i < num_sequences; ++ i) {
    target[i] = new char[strlen(source[i]) + 1];
    strcpy(target[i], source[i]);
  }
  is_sequence_loaded_ = true;
  num_seqs_ = num_sequences;
  return;
}

// reverse the sequences in array "sequence_"
void SequenceBuild::InPlaceReverse(void)  {
  for(int i = 0; i < num_seqs_; ++ i) {
    int l = strlen(sequence_[i]);
    for(int j = 0; j < floor(l / 2); ++ j) {
      swap(sequence_[i][j], sequence_[i][l - 1 - j]);
    }
  }
  return;
}

// building the suffix array on the entire set of sequences (default)
void SequenceBuild::BuildSFADefault(void) {
  this->BuildSuffixArray(sequence_, suffix_array_);
  return;
}

// building the suffix array
void SequenceBuild::BuildSuffixArray(char** sequences, GSA* suffix_array) {
  suffix_array_ = new GSA(sequence_, num_seqs_, true);
  is_sfa_built_ = true;
  return;
}

// building key array that indicates the maximal-extension suffix for each entry
void SequenceBuild::BuildKeyArray(GSA* suffix_array, std::vector<IdxType>& key_array)  {
  IdxType n = suffix_array->getSize();
  key_array.resize(n);
  IdxType block_begin = 0;
  for(IdxType i = 0; i < n - 1; ++ i) {
    // if current suffix length greater than LCP with the next suffix
    // then the current suffix is the end of the block
    if(suffix_array->getSuffixLength(i) > suffix_array->getLcpAt(i + 1)) {
      for(IdxType j = block_begin; j <= i; ++ j) {
        key_array[j] = i;
      }
      block_begin = i + 1;
    }
  }
  // the last block
  for(IdxType j = block_begin; j < n; ++ j) {
    key_array[j] = n - 1;
  }
  return;
}

// building the default key array
void SequenceBuild::BuildKeyArrayDefault(void)  {
  if(!is_sfa_built_)  {
    cout << "Warning: SequenceBuild::BuildKeyArrayDefault attempt to build key array without suffix array built" << endl;
  }
  BuildKeyArray(suffix_array_, key_array_);
  is_k_array_built_ = true;
  return;
} 

// access a given sequence
std::string SequenceBuild::GetSequence(int index)  {
  if(!is_sequence_loaded_)  {
    cout << "Warning: SequenceBuild::GetSequence no sequence loaded" << endl;
    return string("");
  }
  if(index < 0 || index >= num_seqs_)  {
    cout << "Warning: SequenceBuild::GetSequence sequence index out of range" << endl;
    return string("");
  }
  return string(sequence_[index]);
}

std::string SequenceBuild::GetHeader(int index)  {
  if(!is_sequence_loaded_)  {
    cout << "Warning: SequenceBuild::GetSequence no sequence loaded" << endl;
    return string("");
  }
  if(index < 0 || index >= num_seqs_)  {
    cout << "Warning: SequenceBuild::GetSequence sequence index out of range" << endl;
    return string("");
  }
  return string(header_[index]);
}

// access a given suffix from the suffix array
std::string SequenceBuild::GetSuffixSFA(int index) {
  if(!is_sfa_built_)  {
    cout << "Warning: SequenceBuild::GetSuffixSFA suffix array not built" << endl;
    return string("");
  }
  if(index < 0 || index >= suffix_array_->getSize())  {
    cout << "Warning: SequenceBuild::GetSuffixSFA suffix index out of range" << endl;
    return string("");
  }
  return string(suffix_array_->getSuffix(index));
}

// search the suffix array within the object
std::pair<IdxType, IdxType> SequenceBuild::SearchSFA(std::string& search_seed) {
  if(!is_sfa_built_)  {
    cout << "Fatal Error: SequenceBuild::SearchSFA suffix array not built" << endl;
    exit(1);
  }
  return suffix_array_->searchWithLCPs(
      (SfaChar*) search_seed.c_str(), search_seed.length()
  );
}

// find a list of maximal extension reads within the given range
void SequenceBuild::GetMaxExtInfoWithinRange(
    std::pair<IdxType, IdxType>& range, 
    std::list<PositionType>& pos_list
)  {
  if(!is_k_array_built_)  {
    cout << "Fatal Error: SequenceBuild::GetMaxExtRIDWithRange key array not built" << endl;
    exit(1);
  }
  IdxType index = range.first;
  do  {
    index = key_array_[index] + 1;
    PositionType posT;
    posT.rid = (RIDType) suffix_array_->getId(index - 1);
    posT.pos = (POSType) suffix_array_->getPos(index - 1);
    pos_list.push_back(posT);
  } while (index <= range.second);
  return;
}
