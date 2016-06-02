#include "index_sample.h"


using namespace std;

IndexSample::IndexSample(void)  {
  return;
}  

IndexSample::IndexSample(unsigned int in_kmer_size, enum Alphabet in_alphabet)  {
  // copy the kmer size the user-specified
  KmerSize = in_kmer_size;
  // assign the alphabet maps
  init_alphabet(in_alphabet);
  // check parameter setup
  //check_setup();
  KmerHashEntries = 0;
  KmerHashMaxLength = 0;
  return;
}

IndexSample::~IndexSample(void)  {
  return;
}

int IndexSample::AccessAlphabetMap(char c)  {
  if(AlphabetMap.find(c) != AlphabetMap.end())  {
    return AlphabetMap[c];
  }  else  {
    cout << "Warning: The character \'" << c << "\' is not found in the alphabet." << endl;
    return -1;
  }
}

double IndexSample::get_log2_sys_mem(void)  {
  // estimate the available memory
  double log_pages = log2(sysconf(_SC_PHYS_PAGES));
  double log_page_size = log2(sysconf(_SC_PAGE_SIZE));
  return (log_pages + log_page_size);
}

inline void IndexSample::dump_to_HD(const unsigned int cutoff, ostream& idx_fh, ostream& pos_fh)  {
  
  assert(idx_fh.good());
  assert(pos_fh.good());
  
  KmerHashMaxLength = 0;
  vector<unordered_map<KmerType, list<PositionType> >::const_iterator> to_erase;
  for(auto it_hash = KmerHash.begin(); it_hash != KmerHash.end(); ++ it_hash)  {
    if(it_hash->second.size() > cutoff)  {
      for(auto it_vector = it_hash->second.begin(); it_vector != it_hash->second.end(); ++ it_vector)  {
        // write the positions to the position file
        pos_fh.write((char*) &it_vector->rid, sizeof(RIDType));
        pos_fh.write((char*) &it_vector->pos, sizeof(POSType));
        //pos_fh.write((char*) &it_vector, sizeof(PositionType));
      }
      streampos end = pos_fh.tellp();
      // write the start and end position to the indexing file
      idx_fh.write((char *) &it_hash->first, sizeof(KmerType));
      idx_fh.write((char *) &end, sizeof(streampos));
      //cout << (long long) it_hash->first << " " << sizeof(streampos) << " " << (long long) end << endl;
      // collecting memories
      it_hash->second.clear();
      to_erase.resize(to_erase.size() + 1, it_hash);
      -- KmerHashEntries;
    }  else  {
      KmerHashMaxLength = KmerHashMaxLength > it_hash->second.size() ? KmerHashMaxLength : it_hash->second.size();
    }
  }
  
  for(auto it = to_erase.begin(); it != to_erase.end(); ++ it)  {
    KmerHash.erase(*it);
  }
  return;
}

bool IndexSample::ContainSpecialChar(std::string mer)  {
  unsigned int i;
  for(i = 0; i < mer.size(); ++ i) {
    if(AlphabetMap.find(mer[i]) == AlphabetMap.end())  {
      return true;
    } 
  }
  return false;
}

inline void IndexSample::update_position(unsigned int read_id, string seq, ostream& idx_fh, ostream& pos_fh)  {
  
  assert(idx_fh.good());
  assert(pos_fh.good());
  
  unsigned int i;
  for(i = 0; i < seq.length() - KmerSize + 1; ++ i)  {
    string mer = seq.substr(i, KmerSize);
    if(ContainSpecialChar(mer))  {
      // if the k-mer contains special character amino acids that is not contanied in 
      // the 20-alphabet, the k-mer is not used as a seed
      //cout << "special character detected, skiping the kmer \"" << mer << "\" as a seed" << endl;
      continue;
    }
    KmerType encoded_kmer = encode_kmer(mer);
    PositionType position = encode_read_position(read_id, i);
    auto it = KmerHash.find(encoded_kmer);
    if(it != KmerHash.end())  {
      // the kmer already exists
      it->second.push_back(position);
      unsigned s = it->second.size();
      KmerHashMaxLength = KmerHashMaxLength > s ? KmerHashMaxLength : s;
    }  else  {
      // the kmer does not exist, create a new entry
      list<PositionType> list_position;
      list_position.push_back(position);
      KmerHash.insert({{encoded_kmer, list_position}});
      ++ KmerHashEntries;
    }
  }
  // compute the worst case memory consumption: 
  //  64 is the KmerType size, 24 is the vector pointer size, 8 is the PositionType size
  //long int worst_case_mem = KmerHashEntries * (64 + 24 + KmerHashMaxLength * 8);
  // dump it if it is larger than 4G
  //while(worst_case_mem > (long int) pow(2.0, get_log2_sys_mem()))  {
  while(KmerHashEntries > 40000 || KmerHashMaxLength > 10000)  {
    //cout << "dump trigered:  " << KmerHashEntries << " entries, " << KmerHashMaxLength << " records each" << endl;
    // dump all info in the memory into harddist
    //unsigned int cutoff = (int) (KmerHashMaxLength / 2);
    dump_to_HD(0, idx_fh, pos_fh);
    //worst_case_mem = KmerHashEntries * (64 + 24 + KmerHashMaxLength * 8);
  }
  return;
}

void IndexSample::write_format_identifier(string identifier, ostream& file)  {
  // reserve 32 bytes at the beginning of each file
  // the first 30 bytes correspond to the identifier
  // the last 2 bytes correspond to the alphabet in use and the kmer size
   
  assert(identifier.length() < 29);
  assert(file.good() && file.tellp() == 0);
  
  char* out_buff = new char[30];
  for(unsigned int i = 0; i < 30; ++ i)  {
    out_buff[i] = '\0';
  }
  strcpy(out_buff, identifier.c_str());
  char alphabet = (char) CurrentAlphabet;
  char ksize = (char) KmerSize;
  
  file.write((char*) out_buff, 30 * sizeof(char));
  file.write((char*) &alphabet, sizeof(char));
  file.write((char*) &ksize, sizeof(char));
  
  delete [] out_buff;
  
  return;
  
}

// TODO: write read_format_identifier functions
string IndexSample::read_format_identifier(istream& file, enum Alphabet& current_alphabet, unsigned int& kmer_size)  {
  assert(file.good() && file.tellg() == 0);
  char* in_buff = new char[30];
  char ksize_buff, alphabet_buff;
  
  file.read(in_buff, 30 * sizeof(char));
  file.read((char*) &alphabet_buff, sizeof(char));
  file.read((char*) &ksize_buff, sizeof(char));
  
  
  string identifier(in_buff);
  current_alphabet = (enum Alphabet) alphabet_buff;
  kmer_size = (unsigned int) ksize_buff;
  
  delete []in_buff;
  
  return identifier;
}

string IndexSample::read_format_identifier(istream& file)  {
  assert(file.good() && file.tellg() == 0);
  char* in_buff = new char[30];
  char ksize_buff, alphabet_buff;
  
  file.read(in_buff, 30 * sizeof(char));
  file.read((char*) &alphabet_buff, sizeof(char));
  file.read((char*) &ksize_buff, sizeof(char));
  
  string identifier(in_buff);
  
  delete []in_buff;
  
  return identifier;
}

void IndexSample::WriteIndexFile(string in_file, string out_index_file, string out_position_file)  {
  ifstream in_fh(in_file.c_str(), ios_base::in);
  ofstream out_id_fh(out_index_file.c_str(), ios_base::out | ios_base::binary);
  ofstream out_pos_fh(out_position_file.c_str(), ios_base::out | ios_base::binary);
  
  assert(in_fh.good());
  assert(out_id_fh.good());
  assert(out_pos_fh.good());
  
  // write formatting identifier to file
  write_format_identifier("JCVIGUIDEDASSEMBLYINDEX", out_id_fh);
  write_format_identifier("JCVIGUIDEDASSEMBLYPOSITION", out_pos_fh);
  
  unsigned int counter = 0;
  
  string line;
  string appended_seq = "";
  while(getline(in_fh, line))  {
    if(line[0] == '>')  {
      if(appended_seq != "")  {
        update_position(counter - 1, appended_seq, out_id_fh, out_pos_fh);
      }
      ++ counter;
      // clear the concaternated sequence
      appended_seq = "";
    }  else  {
      appended_seq += line;
    }
  }
  // for the last sequence, update the hash
  update_position(counter - 1, appended_seq, out_id_fh, out_pos_fh);
  // enforce a dump for all remaining information
  dump_to_HD(0, out_id_fh, out_pos_fh);
  
  in_fh.close();
  out_id_fh.close();
  out_pos_fh.close();
  
  return;
}

void IndexSample::load_index_file(
  string in_index_file, 
  unordered_map<KmerType, list< pair<streampos, streampos> > >& kmer_positions
)  {
  
  ifstream in_id_fh(in_index_file, ios_base::in | ios_base::binary);
  assert(in_id_fh.good());
  string index_identifier = read_format_identifier(in_id_fh);
  assert(index_identifier == "JCVIGUIDEDASSEMBLYINDEX");
  
  streampos prev = in_id_fh.tellg();
  while(!in_id_fh.eof())  {
    KmerType k_buff;
    streampos e_buff;
    in_id_fh.read((char*) &k_buff, sizeof(KmerType));
    in_id_fh.read((char*) &e_buff, sizeof(streampos));
    kmer_positions[k_buff].push_back({prev, e_buff});
    prev = e_buff;
  }
  in_id_fh.close();
  
  return;
}

inline void IndexSample::load_position_per_kmer(
  istream &pos_fh, 
  unordered_map<KmerType, list< pair<streampos, streampos> > >::const_iterator it,
  list<PositionType>& loaded_positions
)  {
  
  assert(pos_fh.good());
  for(auto it_list = it->second.begin(); it_list != it->second.end(); ++ it_list)  {
    streampos start = it_list->first;
    streampos end = it_list->second;
    assert(end >= start);
    assert((end - start) % (sizeof(RIDType) + sizeof(POSType)) == 0);
    unsigned int num_entries = (int) ((end - start) / (sizeof(RIDType) + sizeof(POSType)));
    // search the right location and start reading
    pos_fh.seekg(start);
    for(unsigned int i = 0; i < num_entries; ++ i)  {
      RIDType rid_buff;
      POSType pos_buff;
      pos_fh.read((char*) &rid_buff, sizeof(RIDType));
      pos_fh.read((char*) &pos_buff, sizeof(POSType));
      PositionType occurrence;
      occurrence.rid = rid_buff;
      occurrence.pos = pos_buff;
      loaded_positions.push_back(occurrence);
    }
  }
  return;
}


void IndexSample::get_position_file_offsets(
  enum Alphabet alphabet,
  unsigned int kmer_size,
  string query_seq, 
  istream& idx_fh, 
  map<streampos, pair<KmerType, streampos> >& file_offsets
)  {
  if(query_seq.length() < kmer_size)  {
    return;
  }
  assert(kmer_size > 0 && kmer_size < 20);
  assert(idx_fh.good() && idx_fh.tellg() == 32);
  
  
  // update the global variable, otherwise will not be able to use many methods
  KmerSize = kmer_size;
  CurrentAlphabet = alphabet;
  init_alphabet(alphabet);
  
  // construct all presented k-mers
  assert(AlphabetMap.size() > 0);
  unsigned int i;
  unordered_map<KmerType, int> presented_kmer_hash;
  for(i = 0; i < query_seq.length() - kmer_size + 1; ++ i)  {
    string s = query_seq.substr(i, kmer_size);
    KmerType key = encode_kmer(s);
    presented_kmer_hash[key] = 1;
  }
  
  // reads the file and retrieve all entries that matches with the presented k-mers
  streampos start = idx_fh.tellg();
  while(!idx_fh.eof())  {
    KmerType k_buff;
    streampos e_buff;
    idx_fh.read((char*) &k_buff, sizeof(KmerType));
    idx_fh.read((char*) &e_buff, sizeof(streampos));
    if(idx_fh.eof())  {
      break;
    }
    // add to offset map if the kmer exists in the presented k-mer hash 
    if(presented_kmer_hash.find(k_buff) != presented_kmer_hash.end())  {
      pair<KmerType, streampos> e(k_buff, e_buff);
      file_offsets[start] = e;
    }
    // update the starting information
    start = e_buff;
  }
  return;
}

void IndexSample::get_position_file_offsets(
  enum Alphabet alphabet,
  unsigned int kmer_size,
  vector<string> query_seq, 
  istream& idx_fh, 
  map<streampos, pair<KmerType, streampos> >& file_offsets
)  {
  assert(kmer_size > 0 && kmer_size < 20);
  assert(idx_fh.good() && idx_fh.tellg() == 32);
  
  
  // update the global variable, otherwise will not be able to use many methods
  KmerSize = kmer_size;
  CurrentAlphabet = alphabet;
  init_alphabet(alphabet);
  
  // construct all presented k-mers
  assert(AlphabetMap.size() > 0);
  unsigned int i, j;
  unordered_map<KmerType, int> presented_kmer_hash;
  for(i = 0; i < query_seq.size(); ++ i) {
    if(query_seq[i].length() < kmer_size)  {
      continue;
    }
    for(j = 0; j < query_seq[i].length() - kmer_size + 1; ++ j)  {
      string s = query_seq[i].substr(j, kmer_size);
      KmerType key = encode_kmer(s);
      presented_kmer_hash[key] = 1;
    }
  }
  // reads the file and retrieve all entries that matches with the presented k-mers
  streampos start = idx_fh.tellg();
  while(!idx_fh.eof())  {
    KmerType k_buff;
    streampos e_buff;
    idx_fh.read((char*) &k_buff, sizeof(KmerType));
    idx_fh.read((char*) &e_buff, sizeof(streampos));
    if(idx_fh.eof())  {
      break;
    }
    // add to offset map if the kmer exists in the presented k-mer hash 
    if(presented_kmer_hash.find(k_buff) != presented_kmer_hash.end())  {
      pair<KmerType, streampos> e(k_buff, e_buff);
      file_offsets[start] = e;
    }
    // update the starting information
    start = e_buff;
  }
  return;
}

void IndexSample::get_kmer_positions(
  istream& pos_fh,
  const map<streampos, pair<KmerType, streampos> >& file_offsets,
  unordered_map<KmerType, list<PositionType> >& kmer_positions
)  {
  assert(pos_fh.good() && pos_fh.tellg() == 32);
  for(auto it = file_offsets.begin(); it != file_offsets.end(); ++ it)  {
    KmerType encoded_kmer = it->second.first;
    streampos start = it->first;
    streampos end = it->second.second;
    assert(end >= start);
    assert((end - start) % (sizeof(RIDType) + sizeof(POSType)) == 0);
    //cout << (long long) encoded_kmer << " " << (long long) start << "  " << (long long) end << endl;
    unsigned int num_entries = (unsigned int) ((end - start) / (sizeof(RIDType) + sizeof(POSType)));
    // search the right location and start reading
    pos_fh.seekg(start);
    for(unsigned int i = 0; i < num_entries; ++ i)  {
      RIDType rid_buff;
      POSType pos_buff;
      pos_fh.read((char*) &rid_buff, sizeof(RIDType));
      pos_fh.read((char*) &pos_buff, sizeof(POSType));
      PositionType occurrence;
      occurrence.rid = rid_buff;
      occurrence.pos = pos_buff;
      //if(kmer_positions.find(encoded_kmer) == kmer_positions.end())  {
      //  
      //}
      kmer_positions[encoded_kmer].push_back(occurrence);
    }
  }
  return;
}

bool IndexSample::IsCompatibleWithSetting(
    const string& idx_file, 
    const string& pos_file
)  {
  ifstream in_id_fh(idx_file, ios_base::in | ios_base::binary);
  ifstream in_pos_fh(pos_file, ios_base::in | ios_base::binary);
  if(!in_id_fh.good() || !in_pos_fh.good())  {
    in_id_fh.close();
    in_pos_fh.close();
    return false;
  }
  enum Alphabet in_use_alphabet_id, in_use_alphabet_pos;
  unsigned int in_use_ksize_id, in_use_ksize_pos;
  string index_identifier = read_format_identifier(in_id_fh, in_use_alphabet_id, in_use_ksize_id);
  string position_identifier = read_format_identifier(in_pos_fh, in_use_alphabet_pos, in_use_ksize_pos);
  if(index_identifier == "JCVIGUIDEDASSEMBLYINDEX" &&
      position_identifier == "JCVIGUIDEDASSEMBLYPOSITION" &&
      in_use_alphabet_id == CurrentAlphabet &&
      in_use_ksize_id == KmerSize &&
      in_use_alphabet_pos == CurrentAlphabet &&
      in_use_ksize_pos == KmerSize)  {
    in_id_fh.close();
    in_pos_fh.close();
    return true;
  }
  in_id_fh.close();
  in_pos_fh.close();
  return false;
}

void IndexSample::ReadKmerIndex(
  string in_index_file, 
  string in_position_file, 
  string query_seq, 
  unordered_map<KmerType, list<PositionType> >& kmer_locations
)  {
  
  ifstream in_id_fh(in_index_file, ios_base::in | ios_base::binary);
  ifstream in_pos_fh(in_position_file, ios_base::in | ios_base::binary);
  assert(in_id_fh.good());
  assert(in_pos_fh.good());
  
  // loads the identifier and parse the in-use alphabet and k-mer size
  enum Alphabet in_use_alphabet_id, in_use_alphabet_pos;
  unsigned int in_use_ksize_id, in_use_ksize_pos;
  
  string index_identifier = read_format_identifier(in_id_fh, in_use_alphabet_id, in_use_ksize_id);
  string position_identifier = read_format_identifier(in_pos_fh, in_use_alphabet_pos, in_use_ksize_pos);
  
  assert(index_identifier == "JCVIGUIDEDASSEMBLYINDEX");
  assert(position_identifier == "JCVIGUIDEDASSEMBLYPOSITION");
  assert(in_use_alphabet_id == in_use_alphabet_pos);
  assert(in_use_ksize_id == in_use_ksize_pos);
  
  // reads in the position information's offsets in the position file
  // map<begin_of_info, pair<encoded_kmer, end_of_info> >
  map<streampos, pair<KmerType, streampos> > offset;   
  get_position_file_offsets(in_use_alphabet_id, in_use_ksize_id, query_seq, in_id_fh, offset);
  // since I have used the sorted map in this case, I will skip the sorting here
  // load read id and positions from the position file
  get_kmer_positions(in_pos_fh, offset, kmer_locations);
  
  in_id_fh.close();
  in_pos_fh.close();
  return;
}

void IndexSample::ReadKmerIndex(
  string in_index_file, 
  string in_position_file, 
  vector<string> query_seq, 
  unordered_map<KmerType, list<PositionType> >& kmer_locations
)  {
  
  ifstream in_id_fh(in_index_file, ios_base::in | ios_base::binary);
  ifstream in_pos_fh(in_position_file, ios_base::in | ios_base::binary);
  assert(in_id_fh.good());
  assert(in_pos_fh.good());
  
  // loads the identifier and parse the in-use alphabet and k-mer size
  enum Alphabet in_use_alphabet_id, in_use_alphabet_pos;
  unsigned int in_use_ksize_id, in_use_ksize_pos;
  
  string index_identifier = read_format_identifier(in_id_fh, in_use_alphabet_id, in_use_ksize_id);
  string position_identifier = read_format_identifier(in_pos_fh, in_use_alphabet_pos, in_use_ksize_pos);
  
  assert(index_identifier == "JCVIGUIDEDASSEMBLYINDEX");
  assert(position_identifier == "JCVIGUIDEDASSEMBLYPOSITION");
  assert(in_use_alphabet_id == in_use_alphabet_pos);
  assert(in_use_ksize_id == in_use_ksize_pos);
  
  // reads in the position information's offsets in the position file
  // map<begin_of_info, pair<encoded_kmer, end_of_info> >
  map<streampos, pair<KmerType, streampos> > offset;   
  get_position_file_offsets(in_use_alphabet_id, in_use_ksize_id, query_seq, in_id_fh, offset);
  // since I have used the sorted map in this case, I will skip the sorting here
  // load read id and positions from the position file
  get_kmer_positions(in_pos_fh, offset, kmer_locations);
  in_id_fh.close();
  in_pos_fh.close();
  return;
}

void IndexSample::BuiltParamLoader(
  const string& in_index_file, const string& in_position_file,
  enum Alphabet& alp_in_use, unsigned int& ksize_in_use
)  {
  ifstream in_id_fh(in_index_file, ios_base::in | ios_base::binary);
  ifstream in_pos_fh(in_position_file, ios_base::in | ios_base::binary);
  assert(in_id_fh.good());
  assert(in_pos_fh.good());
  enum Alphabet in_use_alphabet_id, in_use_alphabet_pos;
  unsigned int in_use_ksize_id, in_use_ksize_pos;
  string index_identifier = read_format_identifier(in_id_fh, in_use_alphabet_id, in_use_ksize_id);
  string position_identifier = read_format_identifier(in_pos_fh, in_use_alphabet_pos, in_use_ksize_pos);
  assert(index_identifier == "JCVIGUIDEDASSEMBLYINDEX");
  assert(position_identifier == "JCVIGUIDEDASSEMBLYPOSITION");
  assert(in_use_alphabet_id == in_use_alphabet_pos);
  assert(in_use_ksize_id == in_use_ksize_pos);
  alp_in_use = in_use_alphabet_id;
  ksize_in_use = in_use_ksize_id;
  in_id_fh.close();
  in_pos_fh.close();
  return;
}

void IndexSample::check_setup()  {
  double expected_log2_mem = log2(sizeof(KmerType)) + KmerSize * log2(AlphabetSize);
  double available_log2_mem = get_log2_sys_mem();
  // if the available mem is not two times larger than the required one, issue a warning
  if(available_log2_mem - expected_log2_mem < 1)  {
    cout << "Memory insufficient. The indexing step may be slow." << endl;
  }  
  return;
}

void IndexSample::init_alphabet(enum Alphabet a)  {
  switch(a)  {
    case ALL20:
      assign_alphabet_ALL20();
      break;
    case DSSP5:
      assign_alphabet_DSSP5();
      break;
    case DSSP10:
      assign_alphabet_DSSP10();
      break;
    case GBMR4:
      assign_alphabet_GBMR4();
      break;
    case GBMR10:
      assign_alphabet_GBMR10();
      break;
    case HSDM5:
      assign_alphabet_HSDM5();
      break;
    case SDM6:
      assign_alphabet_SDM6();
      break;
    case MURPHY5:
      assign_alphabet_MURPHY5();
      break;
    case MURPHY10:
      assign_alphabet_MURPHY10();
      break;
    case TD5:
      assign_alphabet_TD5();
      break;
    case TD10:
      assign_alphabet_TD10();
      break;
    default:
      cout << "The reduced alphabet scheme is not supported!" << endl;
      exit(0);
  }
  return;
}

void IndexSample::assign_alphabet_ALL20()  {
  AlphabetMap = {
    {'P',  0}, {'G',  1}, {'E',  2}, {'K',  3}, {'R',  4}, 
    {'Q',  5}, {'D',  6}, {'S',  7}, {'N',  8}, {'T',  9}, 
    {'H', 10}, {'C', 11}, {'I', 12}, {'V', 13}, {'W', 14}, 
    {'Y', 15}, {'F', 16}, {'A', 17}, {'L', 18}, {'M', 19}
  };
  AlphabetSize = 20;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = ALL20;
  return;
}

void IndexSample::assign_alphabet_DSSP5()  {
  AlphabetMap = {
    {'P',  4}, {'G',  4}, {'E',  0}, {'K',  0}, {'R',  0}, 
    {'Q',  0}, {'D',  3}, {'S',  2}, {'N',  3}, {'T',  2}, 
    {'H',  0}, {'C',  2}, {'I',  1}, {'V',  1}, {'W',  1}, 
    {'Y',  1}, {'F',  1}, {'A',  0}, {'L',  1}, {'M',  1}
  };
  AlphabetSize = 5;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = DSSP5;
  return;
}

void IndexSample::assign_alphabet_DSSP10()  {
  AlphabetMap = {
    {'P',  9}, {'G',  9}, {'E',  0}, {'K',  0}, {'R',  0}, 
    {'Q',  0}, {'D',  8}, {'S',  8}, {'N',  8}, {'T',  6}, 
    {'H',  6}, {'C',  7}, {'I',  1}, {'V',  1}, {'W',  5}, 
    {'Y',  2}, {'F',  3}, {'A',  4}, {'L',  2}, {'M',  4}
  };
  AlphabetSize = 10;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = DSSP10;
  return;
}

void IndexSample::assign_alphabet_GBMR4()  {
  AlphabetMap = {
    {'P',  3}, {'G',  0}, {'E',  1}, {'K',  1}, {'R',  1}, 
    {'Q',  1}, {'D',  1}, {'S',  1}, {'N',  1}, {'T',  1}, 
    {'H',  2}, {'C',  2}, {'I',  2}, {'V',  2}, {'W',  2}, 
    {'Y',  2}, {'F',  2}, {'A',  1}, {'L',  2}, {'M',  2}
  };
  AlphabetSize = 4;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = GBMR4;
  return;
}

void IndexSample::assign_alphabet_GBMR10()  {
  AlphabetMap = {
    {'P',  9}, {'G',  0}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  1}, {'S',  8}, {'N',  2}, {'T',  7}, 
    {'H',  5}, {'C',  6}, {'I',  3}, {'V',  3}, {'W',  3}, 
    {'Y',  4}, {'F',  3}, {'A',  3}, {'L',  3}, {'M',  3}
  };
  AlphabetSize = 10;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = GBMR10;
  return;
}

void IndexSample::assign_alphabet_HSDM5()  {
  AlphabetMap = {
    {'P',  3}, {'G',  3}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  3}, 
    {'H',  4}, {'C',  2}, {'I',  0}, {'V',  0}, {'W',  1}, 
    {'Y',  0}, {'F',  0}, {'A',  3}, {'L',  0}, {'M',  0}
  };
  AlphabetSize = 10;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = HSDM5;
  return;
}

void IndexSample::assign_alphabet_SDM6()  {
  AlphabetMap = {
    {'P',  5}, {'G',  3}, {'E',  3}, {'K',  3}, {'R',  3}, 
    {'Q',  3}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  3}, 
    {'H',  4}, {'C',  1}, {'I',  0}, {'V',  0}, {'W',  2}, 
    {'Y',  0}, {'F',  0}, {'A',  3}, {'L',  0}, {'M',  0}
  };
  AlphabetSize = 6;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = SDM6;
  return;
}

void IndexSample::assign_alphabet_MURPHY5()  {
  AlphabetMap = {
    {'P',  1}, {'G',  1}, {'E',  3}, {'K',  4}, {'R',  4}, 
    {'Q',  3}, {'D',  3}, {'S',  1}, {'N',  3}, {'T',  1}, 
    {'H',  4}, {'C',  0}, {'I',  0}, {'V',  0}, {'W',  2}, 
    {'Y',  2}, {'F',  2}, {'A',  1}, {'L',  0}, {'M',  0}
  };
  AlphabetSize = 5;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = MURPHY5;
  return;
}

void IndexSample::assign_alphabet_MURPHY10()  {
  AlphabetMap = {
    {'P',  8}, {'G',  4}, {'E',  2}, {'K',  1}, {'R',  1}, 
    {'Q',  2}, {'D',  2}, {'S',  9}, {'N',  2}, {'T',  9}, 
    {'H',  5}, {'C',  3}, {'I',  6}, {'V',  6}, {'W',  7}, 
    {'Y',  7}, {'F',  7}, {'A',  0}, {'L',  6}, {'M',  6}
  };
  AlphabetSize = 10;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = MURPHY10;
  return;
}

void IndexSample::assign_alphabet_TD5()  {
  AlphabetMap = {
    {'P',  0}, {'G',  0}, {'E',  1}, {'K',  1}, {'R',  1}, 
    {'Q',  1}, {'D',  2}, {'S',  2}, {'N',  2}, {'T',  2}, 
    {'H',  2}, {'C',  2}, {'I',  3}, {'V',  3}, {'W',  3}, 
    {'Y',  3}, {'F',  3}, {'A',  4}, {'L',  4}, {'M',  4}
  };
  AlphabetSize = 5;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = TD5;
  return;
}

void IndexSample::assign_alphabet_TD10()  {
  AlphabetMap = {
    {'P',  0}, {'G',  1}, {'E',  2}, {'K',  2}, {'R',  2}, 
    {'Q',  2}, {'D',  3}, {'S',  3}, {'N',  3}, {'T',  4}, 
    {'H',  5}, {'C',  5}, {'I',  6}, {'V',  6}, {'W',  7}, 
    {'Y',  7}, {'F',  7}, {'A',  8}, {'L',  9}, {'M',  9}
  };
  AlphabetSize = 10;
  EncodeCharacterBits = ceil(log2(AlphabetSize));
  CurrentAlphabet = TD10;
  return;
}

