#include "../include/sequence_search.h"

using namespace std;

void SequenceSearch::ComputeReducedMap(
    const int mer_len, BioAlphabet &alpha,
    ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
    std::unordered_map<std::string, std::set<KmerUnitType> > &reduced_map
) {
  int i, j;
  KmerUnitcoder coder(alpha, mer_len);
  for(i = 0; i < seqs.size(); ++ i) {
    if(seqs[i].length() < mer_len) continue;
    KmerUnitType ck;
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string re_mer = re_alpha.Convert(seqs[i].substr(j, mer_len));
      if(i == 0)  ck = coder.Encode(seqs[i].substr(0, mer_len).c_str());
      else ck = coder.RightExt(ck, seqs[i][j + mer_len - 1]);
      reduced_map[re_mer].insert(ck);
    }
  }
  return;
}

void SequenceSearch::WriteReducedMap(
    const std::string &file, const int mer_len, BioAlphabet &alpha,
    std::unordered_map<std::string, std::set<KmerUnitType> > &reduced_map
) {
  ofstream out_file;
  out_file.open(file, ios::out);
  if(!out_file.is_open())  {
    cout << "SequenceSearch::DumpSeedMers: Error in writing indexing file " << file << "; Abort." << endl;
  }
  // writing reduced_map
  KmerUnitcoder coder(alpha, mer_len);
  for(auto itr = reduced_map.begin(); itr != reduced_map.end(); ++ itr) {
    out_file << "R:" << itr->first << ":";
    for(auto itrs = itr->second.begin(); itrs != itr->second.end(); ++ itrs) {
      string s = coder.Decode(*itrs);
      out_file << s << ":";
    }
    out_file << endl;
  }
  out_file.close();
  return;
}

void SequenceSearch::IndexReducedMap(
    const int mer_len, BioAlphabet &alpha,
    ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
    const std::string &out_file
)  {
  ofstream out_fh;
  out_fh.open(out_file, ios::out);
  if(!out_fh.is_open())  {
    cout << "SequenceSearch::IndexReducedMap: Error in writing reduced alphabet mapping file " << out_file << "; Abort." << endl;
    exit(1);
  }
  int i, j;
  unordered_map<std::string, std::set<string> > reduced_map;
  int counter = 0;
  for(i = 0; i < seqs.size(); ++ i) {
    if(seqs[i].length() < mer_len) continue;
    KmerUnitType ck;
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mer = seqs[i].substr(j, mer_len);
      string re_mer = re_alpha.Convert(mer);
      reduced_map[re_mer].insert(mer);
      ++ counter;
      if(counter > 1000000)  {
        for(auto itr = reduced_map.begin(); itr != reduced_map.end(); ++ itr) {
          out_fh << "R:" << itr->first << ":";
          for(auto itrs = itr->second.begin(); itrs != itr->second.end(); ++ itrs) {
            out_fh << *itrs << ":";
          }
          out_fh << endl;
        }
        counter = 0;
        reduced_map.clear();
      }
    }
  }
  for(auto itr = reduced_map.begin(); itr != reduced_map.end(); ++ itr) {
    out_fh << "R:" << itr->first << ":";
    for(auto itrs = itr->second.begin(); itrs != itr->second.end(); ++ itrs) {
      out_fh << *itrs << ":";
    }
    out_fh << endl;
  }
  out_fh.close();
  return;
}

void SequenceSearch::LoadReducedMap(
    const std::string &in_file, 
    std::unordered_map<std::string, std::vector<std::string> > &reduced_map
) {
  int i, j;
  string line;
  ifstream in_fh;
  in_fh.open(in_file, ios::in);
  if(!in_fh.is_open())  {
    cout << "SequenceSearch::LoadReducedMap: Error in loading reduced alphabet mapping file " << in_file << "; Abort." << endl;
  }
  // loading the files line by line
  while(getline(in_fh, line))  {
    // extract the information
    vector<int> range;
    range.push_back(0);
    for(i = 0; i < line.length(); ++ i) {
      if(line[i] == ':') range.push_back(i);
    }
    // check if the line corresponds to reduced_alphabet map or sequence ID map
    if(line[0] == 'R')  {
      string key = line.substr(range[1] + 1, range[2] - range[1] - 1);
      for(i = 2; i < range.size() - 1; ++ i) {
        reduced_map[key].push_back(line.substr(range[i] + 1, range[i + 1] - range[i] - 1));
      }
    } else continue;
  }
  in_fh.close();
  return;
}

void SequenceSearch::ComputeSeedMap(
    const int mer_len, ReducedAlphabet &alpha, std::vector<std::string> &seqs,
    std::unordered_map<std::string, std::vector<int> > &seed_mer
) {
  int i, j;
  unordered_map<string, unordered_map<int, bool> > tmp_seed_mer;
  // compute seed_mer
  for(i = 0; i < seqs.size(); ++ i) {
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mer = seqs[i].substr(j, mer_len);
      tmp_seed_mer[mer][i] = true;
    }
  }
  for(auto it = tmp_seed_mer.begin(); it != tmp_seed_mer.end(); ++ it) {
    for(auto iti = it->second.begin(); iti != it->second.end(); ++ iti) {
      seed_mer[it->first].push_back(iti->first);
    }
  }
  return;
}

void SequenceSearch::WriteSeedMap(
    const std::string &file, 
    std::unordered_map<std::string, std::vector<int> > &seed_mer
) {
  ofstream out_file;
  out_file.open(file, ios::app);
  // writing the kmer-readID associations
  for(auto itm = seed_mer.begin(); itm != seed_mer.end(); ++ itm) {
    out_file << "I:" << itm->first << ":";
    for(auto itms = itm->second.begin(); itms != itm->second.end(); ++ itms) {
      out_file << *itms << ":";
    }
    out_file << endl;
  }
  out_file.close();
  return;
}

void SequenceSearch::ComputeKmerFilter(
    const int mer_len, BioAlphabet &alpha, std::vector<std::string> &seqs,
    std::vector<std::string> &all_filters  
) {
  for(auto it = seqs.begin(); it != seqs.end(); ++ it) {
    KmerFiltering kmer_filter(mer_len);
    kmer_filter.KmerBloomFilter(alpha, *it);
    vector<string> ft_str;
    kmer_filter.GenFilterString(ft_str);
    for(auto itf = ft_str.begin(); itf != ft_str.end(); ++ itf) {
      all_filters.push_back(*itf);
    }
  }
  return;
}

void SequenceSearch::WriteKmerFilter(
    const std::string &file, std::vector<std::string> &all_filters
) {
  ofstream out_file;
  out_file.open(file, ios::out);
  // writing the kmer-readID associations
  char c;
  for(auto it = all_filters.begin(); it != all_filters.end(); ++ it) {
    for(int i = 0; i < it->length(); ++ i) {
      c = (*it)[i]; out_file.write((char*) &c, sizeof(char));
    }
  }
  out_file.close();
  return;
}

void SequenceSearch::IndexKmerFilter(
    const int mer_len, BioAlphabet &alpha,
    ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
    const std::string &out_file
) {
  ofstream out_fh;
  out_fh.open(out_file, ios::out);
  if(!out_fh.is_open())  {
    cout << "SequenceSearch::IndexKmerFilter: Error in writing kmer filtering file " << out_file << "; Abort." << endl;
    exit(1);
  }
  vector<string> all_filters;
  for(auto it = seqs.begin(); it != seqs.end(); ++ it) {
    KmerFiltering kmer_filter(mer_len);
    kmer_filter.KmerBloomFilter(alpha, *it);
    vector<string> ft_str;
    kmer_filter.GenFilterString(ft_str);
    kmer_filter.Purge();
    for(auto itf = ft_str.begin(); itf != ft_str.end(); ++ itf) {
      all_filters.push_back(*itf);
    }
    if(all_filters.size() >= 1000) {    
      for(auto it = all_filters.begin(); it != all_filters.end(); ++ it) {
        for(int i = 0; i < it->length(); ++ i) {
          char c = (*it)[i]; out_fh.write((char*) &c, sizeof(char));
        }
      }
      all_filters.clear();
    }  
  }
  for(auto it = all_filters.begin(); it != all_filters.end(); ++ it) {
    for(int i = 0; i < it->length(); ++ i) {
      char c = (*it)[i]; out_fh.write((char*) &c, sizeof(char));
    }
  }
  out_fh.close();
  return;
}

void SequenceSearch::LoadKmerFilter(
    const std::string &file, std::vector<KmerFiltering> &bm_filters
)  {
  vector<string> coded_strings;
  int i, j;
  ifstream in_file;
  in_file.open(file, ios::in);
  if(!in_file.is_open())  {
    cout << "Error: SequenceSearch::LoadKmerFilter: fail to open the indexing file " << file << "; Abort." << endl;
    exit(1);
  }
  // loading the files line by line
  KmerFiltering kmer_param;
  int num_filters = kmer_param.GetNumFilters();
  int filter_size = kmer_param.GetFilterSize();
  int len_filter_str = filter_size / 8; // 8 here is the size of the char type
  string raw_str; char c;
  while(!in_file.eof())  {
    in_file.read((char*) &c, sizeof(char));
    raw_str += c;
  }
  in_file.close();
  int n = raw_str.length() - 1;
  if(n % num_filters == 0 && n % len_filter_str == 0) {
    for(i = 0; i < n; i += len_filter_str) {
      coded_strings.push_back(raw_str.substr(i, len_filter_str));
    }  
  } else  {
    cout << "Error: SequenceSearch::LoadKmerFilter: inconsistent string length, indexing file may be corrupted, please re-build. Abort." << endl;
    exit(1);
  }  
  // construct the bloom filters from the coded strings
  if(coded_strings.size() % num_filters != 0) {
    cout << "Error: SequenceSearch::LoadKmerFilter: inconsistent number of filters, indexing file may be corrupted, please re-build. Abort." << endl;
    exit(1);
  }
  for(i = 0; i < coded_strings.size(); i += num_filters) {
    KmerFiltering kf; kf.LoadFromStringSet(coded_strings, i, i + num_filters - 1);
    bm_filters.push_back(kf);
  }
  return;
}

void SequenceSearch::IndexKmerPosition(
    BioAlphabet &alphabet, const int mer_len, 
    std::vector<std::string> &seqs, const std::string &out_file
) {
  int i, j;
  
  int max_id = pow(alphabet.GetSize(), mer_len) + 1;
  vector<vector<int> > mer_pos;
  mer_pos.resize(max_id + 1);
  KmerUnitcoder coder(alphabet, mer_len);
  
  for(i = 0; i < seqs.size(); ++ i) {
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mr = seqs[i].substr(j, mer_len);
      int mer_id = coder.EncodeInt(mr.c_str());
      mer_pos[mer_id].push_back(i); mer_pos[mer_id].push_back(j); 
    }
  }
  ofstream out_fh;
  out_fh.open(out_file, ios::out);
  // writing the kmer-readID associations
  for(i = 0; i < mer_pos.size(); ++ i) {
    out_fh << i << ":";
    for(j = 0; j < mer_pos[i].size(); ++ j) {
      out_fh << mer_pos[i][j] << ":";
    }
    out_fh << endl;
  }
  out_fh.close();
  return;
}

void SequenceSearch::LoadKmerPosition(
    BioAlphabet &alphabet, const std::string &in_file,
    const int mer_len, std::vector<std::vector<int> >& kmer_pos
) {
  int i, j;
  string line;
  ifstream in_fh;
  in_fh.open(in_file, ios::in);
  if(!in_fh.is_open())  {
    cout << "SequenceSearch::LoadKmerPosition: Error in loading kmer position indexing file: " << in_file << "; Abort." << endl;
    exit(1);
  }
  if(mer_len > 6) {
    cout << "SequenceSearch::LoadKmerPosition: kmer size too large (must be less or equal to 6). Abort." << endl;
    exit(1);
  }
  int max_id = pow(alphabet.GetSize(), mer_len) + 1;
  kmer_pos.resize(max_id + 1);
  KmerUnitcoder coder(alphabet, mer_len);
  // loading the files line by line
  while(getline(in_fh, line))  {
    // extract the information
    vector<int> range;
    range.push_back(0);
    for(i = 0; i < line.length(); ++ i) {
      if(line[i] == ':') range.push_back(i);
    }
    int key_id = stoi(line.substr(range[0], range[1] - range[0]));
    for(i = 1; i < range.size() - 1; ++ i) {
      kmer_pos[key_id].push_back(stoi(line.substr(range[i] + 1, range[i + 1] - range[i] - 1)));
    }
  }
  in_fh.close();
  return;
} 

void SequenceSearch::IndexKmerNeighbor(
    const int mer_len, BioAlphabet &alpha, 
    ScoringProt &fs, const int neighbor_score,
    const std::string &out_file
) {
  int i, j;
  vector<char> all_char;
  for(i = 0; i < alpha.GetSize(); ++ i) {
    all_char.push_back(alpha.GetInvCharMap(i));
  }
  KmerUnitcoder coder(alpha, mer_len);
  // computing all kmers
  vector<string> mer_enum;
  vector<int> mer_id;
  vector<int> cutoff_score;
  mer_enum.resize((int) pow(alpha.GetSize(), mer_len));
  mer_id.resize(mer_enum.size());
  for(i = 0; i < (int) mer_enum.size(); ++ i) {
    string mer(mer_len, ' ');
    int num = i;
    //cout << "num: " << num << endl;
    for(int p = mer_len - 1; p >= 0; -- p) {
      int d = pow(alpha.GetSize(), p);
      int q = (int) (num / d);
      mer[mer_len - p - 1] = all_char[q];
      //cout << " pow, quotion, char:  " << d << " " << q << "  " << alphabet[q] << endl;
      num = num - q * d; 
    }
    mer_enum[i] = mer;
    mer_id[i] = coder.EncodeInt(mer.c_str());
  }  
  // computing all kmer neighbors
  vector<vector<int> > kmer_neighbor;
  kmer_neighbor.resize(mer_enum.size());
  for(i = 0; i < mer_enum.size(); ++ i) {
    kmer_neighbor[mer_id[i]].push_back(mer_id[i]);
    for(j = i + 1; j < mer_enum.size(); ++ j) {
      if(fs.ComputeMatchingScore(mer_enum[i], mer_enum[j]) >= neighbor_score)  {
        kmer_neighbor[mer_id[i]].push_back(mer_id[j]);
        kmer_neighbor[mer_id[j]].push_back(mer_id[i]);
      }
    }
  }
  // writing to output
  ofstream out_fh;
  out_fh.open(out_file, ios::out);
  // writing the kmer-readID associations
  for(i = 0; i < kmer_neighbor.size(); ++ i) {
    // sort the kmers
    sort(kmer_neighbor[i].begin(), kmer_neighbor[i].end());
    out_fh << i << ":";
    for(j = 0; j < kmer_neighbor[i].size(); ++ j) {
      out_fh << kmer_neighbor[i][j] << ":";
    }
    out_fh << endl;
  }
  out_fh.close();
  return;
}

void SequenceSearch::IndexKmerTarget(
    const int mer_len, BioAlphabet &alpha,
    std::vector<std::string> &seqs, 
    ScoringProt &fs, const int neighbor_score,
    const std::string &out_file
)  {
  // construct the kmer neighbor
  int i, j, k;
  vector<char> all_char;
  for(i = 0; i < alpha.GetSize(); ++ i) {
    all_char.push_back(alpha.GetInvCharMap(i));
  }
  KmerUnitcoder coder(alpha, mer_len);
  // computing all kmers
  vector<string> mer_enum;
  vector<int> mer_id;
  vector<int> cutoff_score;
  mer_enum.resize((int) pow(alpha.GetSize(), mer_len));
  mer_id.resize(mer_enum.size());
  for(i = 0; i < (int) mer_enum.size(); ++ i) {
    string mer(mer_len, ' ');
    int num = i;
    //cout << "num: " << num << endl;
    for(int p = mer_len - 1; p >= 0; -- p) {
      int d = pow(alpha.GetSize(), p);
      int q = (int) (num / d);
      mer[mer_len - p - 1] = all_char[q];
      //cout << " pow, quotion, char:  " << d << " " << q << "  " << alphabet[q] << endl;
      num = num - q * d; 
    }
    mer_enum[i] = mer;
    mer_id[i] = coder.EncodeInt(mer.c_str());
  }  
  // computing all kmer neighbors
  vector<vector<int> > kmer_neighbor;
  kmer_neighbor.resize(mer_enum.size());
  for(i = 0; i < mer_enum.size(); ++ i) {
    for(j = i; j < mer_enum.size(); ++ j) {
      if(fs.ComputeMatchingScore(mer_enum[i], mer_enum[j]) >= neighbor_score)  {
        kmer_neighbor[mer_id[i]].push_back(mer_id[j]);
        kmer_neighbor[mer_id[j]].push_back(mer_id[i]);
      }
    }
  }
  // estimate the target region size
  vector<int> target_count(kmer_neighbor.size(), 0);
  for(i = 0; i < seqs.size(); ++ i) {
    for(j = 0; j < kmer_neighbor.size(); ++ j) {
      ++ target_count[j];
    }
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mr = seqs[i].substr(j, mer_len);
      int mer_id = coder.EncodeInt(mr.c_str());
      for(k = 0; k < kmer_neighbor[mer_id].size(); ++ k) {
        ++ target_count[kmer_neighbor[mer_id][k]];
      }
    }
  }
  // store the information
  short **kmer_target = new short* [kmer_neighbor.size()];
  for(i = 0; i < kmer_neighbor.size(); ++ i) {
    kmer_target[i] = new short [target_count[i] + 1];
  }
  vector<int> target_index(kmer_neighbor.size(), 0);
  for(i = 0; i < seqs.size(); ++ i) {
    for(j = 0; j < kmer_neighbor.size(); ++ j) {
      kmer_target[j][target_index[j]] = (short) -1;
      ++ target_index[j];
    }
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mr = seqs[i].substr(j, mer_len);
      int mer_id = coder.EncodeInt(mr.c_str());
      for(k = 0; k < kmer_neighbor[mer_id].size(); ++ k) {
        int id = kmer_neighbor[mer_id][k];
        kmer_target[id][target_index[id]] = (short) j;
        ++ target_index[id];
      }
    }
  }
  // collect memory
  for(i = 0; i < kmer_neighbor.size(); ++ i) {
    delete [] kmer_target[i];
  }
  delete [] kmer_target;
  return;
}

void SequenceSearch::LoadKmerNeighbor(
    BioAlphabet &alphabet, const std::string &in_file, const int mer_len,
    std::vector<std::vector<int> >& kmer_neighbor
) {
  int i, j;
  string line;
  ifstream in_fh;
  in_fh.open(in_file, ios::in);
  if(!in_fh.is_open())  {
    cout << "SequenceSearch::LoadKmerPosition: Error in loading kmer neighbhor indexing file: " << in_file << "; Abort." << endl;
    exit(1);
  }
  int max_id = pow(alphabet.GetSize(), mer_len) + 1;
  kmer_neighbor.resize(max_id + 1);
  KmerUnitcoder coder(alphabet, mer_len);
  // loading the files line by line
  while(getline(in_fh, line))  {
    // extract the information
    vector<int> range;
    range.push_back(0);
    for(i = 0; i < line.length(); ++ i) {
      if(line[i] == ':') range.push_back(i);
    }
    int key_id = stoi(line.substr(range[0], range[1] - range[0]));
    for(i = 1; i < range.size() - 1; ++ i) {
      int neighbor_id = stoi(line.substr(range[i] + 1, range[i + 1] - range[i] - 1));
      kmer_neighbor[key_id].push_back(neighbor_id);
    }
  }
  in_fh.close();
  return;
} 

void SequenceSearch::ExtractSeedMers(
    const int mer_len, ReducedAlphabet &alpha, std::vector<std::string> &seqs,
    std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
    std::unordered_map<std::string, std::vector<int> > &seed_mer
) {
  int i, j;
  
  unordered_map<string, unordered_map<int, bool> > tmp_seed_mer;
  // compute seed_mer
  for(i = 0; i < seqs.size(); ++ i) {
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mer = seqs[i].substr(j, mer_len);
      tmp_seed_mer[mer][i] = true;
    }
  }
  for(auto it = tmp_seed_mer.begin(); it != tmp_seed_mer.end(); ++ it) {
    for(auto iti = it->second.begin(); iti != it->second.end(); ++ iti) {
      seed_mer[it->first].push_back(iti->first);
    }
  }
  // compute reduced_map
  for(auto it = seed_mer.begin(); it != seed_mer.end(); ++ it) {
    reduced_map[alpha.Convert(it->first)].push_back(it->first);
  }
  return;
}

void SequenceSearch::DumpSeedMers(
    const std::string &file,
    std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
    std::unordered_map<std::string, std::vector<int> > &seed_mer
) {
  ofstream out_file;
  out_file.open(file, ios::out);
  if(!out_file.is_open())  {
    cout << "SequenceSearch::DumpSeedMers: Error in writing indexing file " << file << "; Abort." << endl;
  }
  // writing reduced_map
  for(auto itr = reduced_map.begin(); itr != reduced_map.end(); ++ itr) {
    out_file << "R:" << itr->first << ":";
    for(auto itrs = itr->second.begin(); itrs != itr->second.end(); ++ itrs) {
      out_file << *itrs << ":";
    }
    out_file << endl;
  }
  // writing the kmer-readID associations
  for(auto itm = seed_mer.begin(); itm != seed_mer.end(); ++ itm) {
    out_file << "I:" << itm->first << ":";
    for(auto itms = itm->second.begin(); itms != itm->second.end(); ++ itms) {
      out_file << *itms << ":";
    }
    out_file << endl;
  }
  out_file.close();
  return;
}

void SequenceSearch::LoadSeedMers(
    const std::string &file,
    std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
    std::unordered_map<std::string, std::vector<int> > &seed_mer
) {
  int i, j;
  string line;
  ifstream in_file;
  in_file.open(file, ios::in);
  if(!in_file.is_open())  {
    cout << "SequenceSearch::LoadSeedMers: Error in loading indexing file " << file << "; Abort." << endl;
  }
  // loading the files line by line
  while(getline(in_file, line))  {
    // extract the information
    vector<int> range;
    range.push_back(0);
    for(i = 0; i < line.length(); ++ i) {
      if(line[i] == ':') range.push_back(i);
    }
    // check if the line corresponds to reduced_alphabet map or sequence ID map
    if(line[0] == 'R')  {
      string key = line.substr(range[1] + 1, range[2] - range[1] - 1);
      for(i = 2; i < range.size() - 1; ++ i) {
        reduced_map[key].push_back(line.substr(range[i] + 1, range[i + 1] - range[i] - 1));
      }
    } else if(line[0] == 'I') {
      string key = line.substr(range[1] + 1, range[2] - range[1] - 1);
      for(i = 2; i < range.size() - 1; ++ i) {
        seed_mer[key].push_back(stoi(line.substr(range[i] + 1, range[i + 1] - range[i] - 1)));
      }
    }
  }
  in_file.close();
  
  //*****************************
  //for(auto it = reduced_map.begin(); it != reduced_map.end(); ++ it) {
  //  cout << "***: " << it->second.size() << endl;
  //}
  //*****************************
  
  return;
}

bool SortAlignInterval(const AlignIntervalType &i1, const AlignIntervalType &i2)  {
  if(i1.q1 < i2.q1 || (i1.q1 == i2.q1 && i1.t1 < i2.t1))  {
    return true;
  }
  return false;
}

void SequenceSearch::SelectID(
    BioAlphabet &alphabet, const int mer_len, 
    const int screen_score, std::string &query, 
    std::vector<std::string> &seqs, ScoringProt &fs,
    std::vector<std::vector<int> >& kmer_neighbor,
    std::vector<std::pair <int, int> > &q_interval, 
    std::vector<std::pair <int, int> > &t_interval,
    std::vector<int> &mx_score, std::vector<int> &t_ID
) {
  // compute high-score matches for each kmer in the group with the same reduce alphabet string
  int i, j, k, l;
  int max_id = pow(alphabet.GetSize(), mer_len) + 1;
  KmerUnitcoder coder(alphabet, mer_len); 
  int mer_id; 
  
  int max_len = 0;
  for(i = 0; i < seqs.size(); ++ i) {
    max_len = max_len > seqs[i].length() ? max_len : seqs[i].length();
  }
  
  vector<vector <int> > mer_q_loc;
  mer_q_loc.resize(max_id);
  for(i = 0; i <= query.length() - mer_len; ++ i) {
    if(i == 0)  mer_id = coder.EncodeInt(query.substr(0, mer_len).c_str());
    else  mer_id = coder.EncodeIntRight(mer_id, query[i + mer_len - 1]);
    for(j = 0; j < kmer_neighbor[mer_id].size(); ++ j) {
      mer_q_loc[kmer_neighbor[mer_id][j]].push_back(i);
    }
  }
  
  int asize = 2 * (query.length() + max_len), abit = asize * sizeof(int);
  int *diff_array = new int [asize];
  for(i = 0; i < seqs.size(); ++ i) {
    //diff_array = vector<int>(2 * (query.length() + max_len), -1);
    memset(diff_array, -1, abit);
    int mx_hsp_score = 0;
    int mx_q1, mx_q2, mx_t1, mx_t2;
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      if(j == 0) mer_id = coder.EncodeInt(seqs[i].substr(0, mer_len).c_str());
      else mer_id = coder.EncodeIntRight(mer_id, seqs[i][j + mer_len - 1]);
      for(k = 0; k < mer_q_loc[mer_id].size(); ++ k) {
      
        int diff = j - mer_q_loc[mer_id][k] + query.length();
        int idx = diff * 2;
        if(diff_array[idx] == -1)  {
          diff_array[idx] = mer_q_loc[mer_id][k]; diff_array[idx + 1] = j;
        } else {
          if(mer_q_loc[mer_id][k] - diff_array[idx] == j - diff_array[idx + 1] && 
             j - diff_array[idx + 1] > mer_len) {
            
            int q1l = diff_array[idx], q2r = mer_q_loc[mer_id][k] + mer_len - 1;
            int t1l = diff_array[idx + 1], t2r = j + mer_len - 1;            
            int hsp_score = fs.ComputeMatchingScoreExtend(
                query, seqs[i], 25, q1l, q2r, t1l, t2r
            );           
            if(hsp_score >= mx_hsp_score) {
              mx_hsp_score = hsp_score;
              mx_q1 = q1l; mx_q2 = q2r; 
              mx_t1 = t1l; mx_t2 = t2r;  
            }
            // update the array
            diff_array[idx] = mer_q_loc[mer_id][k]; diff_array[idx + 1] = j;
          }
        }
      }
    }
    if(mx_hsp_score >= screen_score)  {
      mx_score.push_back(mx_hsp_score);
      t_ID.push_back(i);
      q_interval.push_back(make_pair(mx_q1, mx_q2));
      t_interval.push_back(make_pair(mx_t1, mx_t2));
      //cout << "score: " << mx_hsp_score << " ~ ID: " << i << "  " << mx_q1 << "  " << mx_q2 << " " << mx_t1 << " " << mx_t2 << endl;
    }
  }
  delete [] diff_array;
  return;
}


int SequenceSearch::SelectID(
    BioAlphabet &alphabet, const int mer_len, 
    const int screen_score, std::string &query, 
    std::vector<std::string> &seqs, ScoringProt &fs,
    std::vector<std::vector<int> >& kmer_neighbor,
    int *q_interval, int *t_interval,
    int *mx_score, int *t_ID
) {
  // compute high-score matches for each kmer in the group with the same reduce alphabet string
  int i, j, k, l;
  int max_id = pow(alphabet.GetSize(), mer_len) + 1;
  KmerUnitcoder coder(alphabet, mer_len); 
  int mer_id; 
  
  int max_len = 0;
  for(i = 0; i < seqs.size(); ++ i) {
    max_len = max_len > seqs[i].length() ? max_len : seqs[i].length();
  }
  
  vector<vector <int> > mer_q_loc;
  mer_q_loc.resize(max_id);
  for(i = 0; i <= query.length() - mer_len; ++ i) {
    if(i == 0)  mer_id = coder.EncodeInt(query.substr(0, mer_len).c_str());
    else  mer_id = coder.EncodeIntRight(mer_id, query[i + mer_len - 1]);
    for(j = 0; j < kmer_neighbor[mer_id].size(); ++ j) {
      mer_q_loc[kmer_neighbor[mer_id][j]].push_back(i);
    }
  }
  //cout << "Done using kmer-neighbor" << endl;
  
  int num_cand = 0;
  int asize = 2 * (query.length() + max_len), abit = asize * sizeof(int);
  int *diff_array = new int [asize];
  for(i = 0; i < seqs.size(); ++ i) {
    //diff_array = vector<int>(2 * (query.length() + max_len), -1);
    memset(diff_array, -1, abit);
    int mx_hsp_score = 0;
    int mx_q1, mx_q2, mx_t1, mx_t2;
    char *t_seq = new char [seqs[i].length() + 1];
    strcpy(t_seq, seqs[i].c_str());
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      
      //if(j == 0) mer_id = coder.EncodeInt(seqs[i].substr(0, mer_len).c_str());
      //else mer_id = coder.EncodeIntRight(mer_id, seqs[i][j + mer_len - 1]);
      
      if(j == 0) mer_id = coder.EncodeInt(t_seq);
      else mer_id = coder.EncodeIntRight(mer_id, t_seq[j + mer_len - 1]);
      
      
      for(k = 0; k < mer_q_loc[mer_id].size(); ++ k) {
      
        int diff = j - mer_q_loc[mer_id][k] + query.length();
        int idx = diff * 2;
        if(diff_array[idx] == -1)  {
          diff_array[idx] = mer_q_loc[mer_id][k]; diff_array[idx + 1] = j;
        } else {
          if(mer_q_loc[mer_id][k] - diff_array[idx] == j - diff_array[idx + 1] && 
             j - diff_array[idx + 1] > mer_len) {
            
            int q1l = diff_array[idx], q2r = mer_q_loc[mer_id][k] + mer_len - 1;
            int t1l = diff_array[idx + 1], t2r = j + mer_len - 1;            
            int hsp_score = fs.ComputeMatchingScoreExtend(
                query, seqs[i], 25, q1l, q2r, t1l, t2r
            );           
            if(hsp_score >= mx_hsp_score) {
              mx_hsp_score = hsp_score;
              mx_q1 = q1l; mx_q2 = q2r; 
              mx_t1 = t1l; mx_t2 = t2r;  
            }
            // update the array
            diff_array[idx] = mer_q_loc[mer_id][k]; diff_array[idx + 1] = j;
          }
        }
      }
      
    }
    delete [] t_seq;
    //cout << "Done filtering" << endl;
    if(mx_hsp_score >= screen_score)  {
      mx_score[num_cand] = mx_hsp_score;
      t_ID[num_cand] = i;
      q_interval[2 * num_cand] = mx_q1;
      q_interval[2 * num_cand + 1] = mx_q2;
      t_interval[2 * num_cand] = mx_t1;
      t_interval[2 * num_cand + 1] = mx_t2;
      ++ num_cand;
    }
    //cout << "Done writing information" << endl;
  }
  delete [] diff_array;
  
  return num_cand;
}

void SequenceSearch::SelectID(
    const int mer_len, BioAlphabet &alpha, ReducedAlphabet &re_alpha, 
    std::string &query, std::vector<std::string> &seqs,
    std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
    std::vector<KmerFiltering> &bm_filters, double scale, ScoringProt &fs, 
    std::vector<std::vector<int> > &seed_positions
) {
  // compute high-score matches for each kmer in the group with the same reduce alphabet string
  int i, j, k;
  unordered_map<string, int> q_mer;
  for(i = 0; i <= query.length() - mer_len; ++ i) {
    string mer = query.substr(i, mer_len);
    string re_mer = re_alpha.Convert(mer);
    int cutoff = fs.GetMaxScore(mer) * scale;
    if(reduced_map.find(re_mer) == reduced_map.end()) continue;
    for(auto it = reduced_map[re_mer].begin(); it != reduced_map[re_mer].end(); ++ it) {
      if(fs.ComputeMatchingScore(mer, *it) >= cutoff)  {
        q_mer[*it] = i;
      }
    }
  }
  // computing the unitigs that contains kmers being hashed for each position
  KmerFiltering filter_check(mer_len);
  vector<vector<int> > q_mer_filters; vector<string> q_mer_array;
  for(auto it = q_mer.begin(); it != q_mer.end(); ++ it) {
    vector<int> check_pos;
    filter_check.ComputeSingleKmerHash(alpha, it->first, check_pos);
    q_mer_filters.push_back(check_pos); q_mer_array.push_back(it->first);
  }
  // checking whether the kmer is present in a sequence 
  bool *taken_seqs = new bool [bm_filters.size()];
  memset(taken_seqs, false, bm_filters.size());
  for(i = 0; i < bm_filters.size(); ++ i) {
    for(j = 0; j < q_mer_filters.size(); ++ j) {
      //cout << "Check bloom filter:  " << q_mer_array[j] << "  " << bm_filters[i].IsKmerUnitPresent(q_mer_filters[j]) << " " << bm_filters[i].IsKmerUnitPresent(coder_debug.Encode(q_mer_array[j].c_str())) << " " << checker_debug.IsKmerUnitPresent(coder_debug.Encode(q_mer_array[j].c_str())) << endl;
      taken_seqs[i] = taken_seqs[i] | bm_filters[i].IsKmerUnitPresent(q_mer_filters[j]);
    }
  }
  // compute seeding information
  seed_positions.resize(bm_filters.size());
  for(i = 0; i < bm_filters.size(); ++ i) {
    if(!taken_seqs[i])  {continue;}       
    for(j = 0; j <= seqs[i].length() - mer_len; ++ j) {
      string mer = seqs[i].substr(j, mer_len);
      if(q_mer.find(mer) != q_mer.end())  {
        seed_positions[i].push_back(q_mer[mer]);
        seed_positions[i].push_back(j);
      }
    }
  }
  delete [] taken_seqs;
  return;
}


void SequenceSearch::DefineChainingPreSet(
    const int band_size,
    const std::vector<AlignIntervalType> &matching, 
    std::vector<std::vector<int> > &pre_set
) {
  int i, j;
  pre_set.resize(matching.size());
  for(i = 0; i < matching.size() - 1; ++ i) {
    int q_bound = MAX, t_bound = MAX;
    for(j = i + 1; j < matching.size(); ++ j) {
      if(matching[j].q1 <= matching[i].q2 || matching[j].t1 <= matching[i].t2)  continue;
      int diff_q = matching[j].q1 - matching[i].q1;
      int diff_t = matching[j].t1 - matching[i].t1;
      if((matching[j].q1 < q_bound || matching[j].t1 < t_bound) && 
         (abs(diff_q - diff_t) <= band_size)
      )  {
        q_bound = q_bound < matching[j].q2 ? q_bound : matching[j].q2;
        t_bound = t_bound < matching[j].t2 ? t_bound : matching[j].t2;
        pre_set[j].push_back(i);
      }
      if(matching[j].q1 > q_bound && matching[j].t1 > t_bound)  break;
    }
  }
  
  /*
  for(i = 0; i < matching.size(); ++ i) {
    cout << "(" << matching[i].q1 << "  " << matching[i].q2 << "  " << matching[i].t1 << "  " << matching[i].t2 << ")" << endl;
  }
  
  for(i = 0; i < pre_set.size(); ++ i) {
    for(j = 0; j < pre_set[i].size(); ++ j) {
      cout << pre_set[i][j] << "  ";
    }
    cout << endl;
  }
  */
  return;
}

void SequenceSearch::ExtendSeedToMatching(
    const int mer_len, const std::string& query, const std::string &target, 
    const std::vector<int> &all_seeds, ScoringProt &fs,
    std::vector<AlignIntervalType> &matching 
) {
  // record the intervals
  vector<AlignIntervalType> seed_matching; 
  for(int i = 0; i < all_seeds.size(); i += 2) {
    AlignIntervalType x; 
    x.q1 = all_seeds[i]; x.q2 = x.q1 + mer_len - 1;
    x.t1 = all_seeds[i + 1]; x.t2 = x.t1 + mer_len - 1;
    //x.score = fs.ComputeMatchingScore(query, x.q1, target, x.t1, mer_len);
    //x.score = fs.ComputeMatchingScoreExtend(query, target, x.q1, x.q2, x.t1, x.t2);
    //cout << "Pair score:  " << x.score << endl;
    seed_matching.push_back(x);
  }
  //sort(seed_matching.begin(), seed_matching.end(), SortAlignInterval);
  
  /*
  for(int i = 0; i < seed_matching.size(); ++ i) {
    cout << seed_matching[i].q1 << "  " << seed_matching[i].t1 << " ";
  }
  cout << endl;
  */
  // extend the intervals one by one
  vector<bool> valid_matching(seed_matching.size(), true);
  for(int i = 0; i < seed_matching.size(); ++ i) {
    if(!valid_matching[i]) continue;
    seed_matching[i].score = fs.ComputeMatchingScoreExtend(
        query, target, 50, 
        seed_matching[i].q1, seed_matching[i].q2,
        seed_matching[i].t1, seed_matching[i].t2
    );
    for(int j = i + 1; j < seed_matching.size(); ++ j) {
      if(seed_matching[j].q2 <= seed_matching[i].q2 &&
         seed_matching[j].t2 <= seed_matching[i].t2 &&
         seed_matching[i].q2 - seed_matching[j].q2 == seed_matching[i].t2 - seed_matching[j].t2)  {
        valid_matching[j] = false;
      }
      if(seed_matching[j].q1 > seed_matching[i].q2) break;
    }
  }
  // copy the valid matchings
  for(int i = 0; i < seed_matching.size(); ++ i) {
    if(valid_matching[i]) matching.push_back(seed_matching[i]);
  }
  return;
}

void SequenceSearch::ExtendSeedToMatchingTwoHit(
    const int mer_len, const std::string& query, const std::string &target, 
    const std::vector<int> &all_seeds, ScoringProt &fs,
    std::vector<AlignIntervalType> &matching 
) {
  // record the intervals
  vector<AlignIntervalType> seed_matching; 
  for(int i = 0; i < all_seeds.size(); i += 4) {
    AlignIntervalType x; 
    x.q1 = all_seeds[i]; x.q2 = all_seeds[i + 2] + mer_len - 1;
    x.t1 = all_seeds[i + 1]; x.t2 = all_seeds[i + 3] + mer_len - 1;
    //cout << "LLL: " << query.length() << "  " << target.length() << endl; 
    //cout << "AAA: " << all_seeds[i] << "  " << all_seeds[i + 1] << "  " << all_seeds[i + 2] << "  " << all_seeds[i + 3] << endl;
    //cout << "BBB: " << x.q1 << "  " << x.q2 << "  " << x.t1 << "  " << x.t2 << endl;
    //x.score = fs.ComputeMatchingScore(query, x.q1, target, x.t1, mer_len);
    //x.score = fs.ComputeMatchingScoreExtend(query, target, x.q1, x.q2, x.t1, x.t2);
    //cout << "Pair score:  " << x.score << endl;
    seed_matching.push_back(x);
  }
  //sort(seed_matching.begin(), seed_matching.end(), SortAlignInterval);
  
  /*
  for(int i = 0; i < seed_matching.size(); ++ i) {
    cout << seed_matching[i].q1 << "  " << seed_matching[i].t1 << " ";
  }
  cout << endl;
  */
  // extend the intervals one by one
  vector<bool> valid_matching(seed_matching.size(), true);
  for(int i = 0; i < seed_matching.size(); ++ i) {
    if(!valid_matching[i]) continue;
    seed_matching[i].score = fs.ComputeMatchingScoreExtend(
        query, target, 25, 
        seed_matching[i].q1, seed_matching[i].q2,
        seed_matching[i].t1, seed_matching[i].t2
    );
    for(int j = i + 1; j < seed_matching.size(); ++ j) {
      if(seed_matching[j].q2 <= seed_matching[i].q2 &&
         seed_matching[j].t2 <= seed_matching[i].t2 &&
         seed_matching[i].q2 - seed_matching[j].q2 == seed_matching[i].t2 - seed_matching[j].t2)  {
        valid_matching[j] = false;
      }
      if(seed_matching[j].q1 > seed_matching[i].q2) break;
    }
  }
  // copy the valid matchings
  for(int i = 0; i < seed_matching.size(); ++ i) {
    if(valid_matching[i]) matching.push_back(seed_matching[i]);
  }
  return;
}

int SequenceSearch::FindBestMatching(
    const int mer_len, const int screen_score,
    const std::string& query, const std::string &target, 
    const std::vector<AlignIntervalType> &matching, ScoringProt &fs, 
    std::vector<std::pair<int, int> > &q_interval, 
    std::vector<std::pair<int, int> > &t_interval,
    std::vector<int> &mx_pair_score
) {
  if(matching.size() <= 0)  {
    cout << "Warning: SequenceSearch::FindBestMatching: empty input possible matching set!" << endl;
    pair<int, int> foo; foo.first = foo.second = -1;
    q_interval.push_back(foo); t_interval.push_back(foo);
    int ws = 2 * fs.gap_open_ + fs.gap_ext_ * (query.length() + target.length());
    mx_pair_score.push_back(ws);
    return ws;  
  }
  int i, j, k;
  
  vector<vector<int> > pre_set;
  DefineChainingPreSet(20, matching, pre_set);
  AlignBatch aligner;
  // using DP to find the set of best matching
  vector<int> max_score(matching.size(), 0);
  vector<int> prev_intv(matching.size(), -1);
  vector<bool> taken(matching.size(), false);  
  // using DP to compute the optimal non-overlapping set for the query
  for(i = 0; i < matching.size(); ++ i) {
    max_score[i] = matching[i].score; 
    prev_intv[i] = -1; taken[i] = true;
    for(j = 0; j < pre_set[i].size(); ++ j) {
      k = pre_set[i][j];
      if(matching[k].q2 < matching[i].q1 && matching[k].t2 < matching[i].t1)  {
        string gs = query.substr(matching[k].q2 + 1, matching[i].q1 - matching[k].q2 - 1);
        string gt = target.substr(matching[k].t2 + 1, matching[i].t1 - matching[k].t2 - 1);
        int gb = abs((int) gs.length() - (int) gt.length()) * 2 + 10;
        int ns = max_score[k] + matching[i].score
            + aligner.AlignGlobalPairwise(gs, gt, gb, fs);
            //+ aligner.AlignGlobalSingleEditDist(gs, gt, gb);
        if(ns > max_score[i])  {
          max_score[i] = ns; prev_intv[i] = k; taken[i] = true;
        } 
      }
    }
  }
  // trace-back and identify the set of sequences to be aligned
  int num_matchings = 0, current_idx, tail_idx, head_idx;
  unordered_map<int, bool> taken_pair;
  
  int mx;
  //do {
  mx = 0;
  for(i = 0; i < max_score.size(); ++ i) {
    if(taken_pair.find(i) == taken_pair.end() && max_score[i] > mx) {
      mx = max_score[i]; current_idx = i;
    }
  }
  //cout << "max location:  " << current_idx << endl;
  //if(mx <= 0 || taken_pair.find(current_idx) != taken_pair.end())  { break;  }
  tail_idx = current_idx;
  while(current_idx >= 0) {
    if(taken[current_idx])  {
      ++ num_matchings; head_idx = current_idx; taken_pair[current_idx] = true;
      //cout << "taken: " << current_idx << endl;
    }
    current_idx = prev_intv[current_idx];
  }
  pair<int, int> qv, tv;
  qv.first = matching[head_idx].q1; qv.second = matching[tail_idx].q2;
  tv.first = matching[head_idx].t1; tv.second = matching[tail_idx].t2;
  //cout << "Num pairs taken: " << np << endl;
  //if(mx >= screen_score)  {
  q_interval.push_back(qv);
  t_interval.push_back(tv);
  mx_pair_score.push_back(mx);
  //}
  //cout << "Final: " << qv.first << "  " << qv.second << " " << tv.first << "  " << tv.second << endl;
  //} while(mx >= screen_score);
  return num_matchings;
}

void SequenceSearch::SelectAlignmentRegions(
    const int mer_len, const int screen_score, ScoringProt &fs,
    std::string &query, std::vector<std::string> &seqs, 
    std::vector<std::vector<int> > &seed_positions,
    std::vector<std::pair <int, int> > &q_interval, 
    std::vector<std::pair <int, int> > &t_interval,
    std::vector<int> &mx_score, std::vector<int> &t_ID    
) {
  int m, i, j; 
  if(seed_positions.size() != seqs.size())  {
    cout << "Error: SequenceSearch::SelectAlignmentRegions: Unmatched size of seeding information and target sequences, abort." << endl;
    exit(1);
  }
  // for each target sequence compute the corresponding alignment score
  for(m = 0; m < seed_positions.size(); ++ m) {
    // attemp to filter target sequences based on kmer location information
    if(seed_positions[m].size() <= 0) continue;
    
    /*
    cout << "***********************" << endl;
    for(i = 0; i < seed_positions[m].size(); ++ i)  {
      cout << seed_positions[m][i] << " ";
    }
    cout << endl;
    */
    
    /*
    vector<int> seed_pair;
    bool is_passed = false;
    for(i = 0; i < seed_positions[m].size() - 2; i += 2) {
      for(j = i + 2; j < seed_positions[m].size(); j += 2)  {
        if(seed_positions[m][j] - seed_positions[m][i] > 20) break;
        int x = seed_positions[m][j] - seed_positions[m][i];
        int y = seed_positions[m][j + 1] - seed_positions[m][i + 1];
        if(x >= mer_len && x == y)  {
          //is_passed = true;
          //break;
          seed_pair.push_back(seed_positions[m][i]);
          seed_pair.push_back(seed_positions[m][i + 1]);
          seed_pair.push_back(seed_positions[m][j]);
          seed_pair.push_back(seed_positions[m][j + 1]);
        }
      }
      //if(is_passed) break;
    }
    if(seed_pair.size() <= 0) {
      //cout << "filtered" << endl;
      continue;
    }
    */
    
    // compute rough alignment between the query and the target
    //cout << "Sequence ID: " << m << endl;
    vector<AlignIntervalType> matching;
    //ExtendSeedToMatching(mer_len, query, seqs[m], seed_positions[m], fs, matching);
    //cout << "Before calling ExtendTwoHit" << endl;
    //ExtendSeedToMatchingTwoHit(mer_len, query, seqs[m], seed_pair, fs, matching);
    ExtendSeedToMatchingTwoHit(mer_len, query, seqs[m], seed_positions[m], fs, matching);
    //cout << "End calling ExtendTwoHit" << endl;
    //cout << "size:  " << matching.size() << endl;
    int max_matching_score = 0;
    int max_index = 0;
    for(i = 0; i < matching.size(); ++ i) {
      if(matching[i].score > max_matching_score)  {
        max_matching_score = matching[i].score;
        max_index = i;
      }
    }
    
    int q1 = matching[max_index].q1, q2 = matching[max_index].q2;
    int t1 = matching[max_index].t1, t2 = matching[max_index].t2;   
    int pre_len = q1 < t1 ? q1 : t1;
    int post_len = query.length() - q2 < seqs[m].length() - t2 ?
        query.length() - q2 : seqs[m].length() - t2;
    q_interval.push_back(make_pair(q1 - pre_len, q2 + post_len - 1));
    t_interval.push_back(make_pair(t1 - pre_len, t2 + post_len - 1));
    t_ID.push_back(m);
    mx_score.push_back(max_matching_score);
  }
  return;
}

void SequenceSearch::FetchAlignmentSeqs(
    std::string &query, std::vector<std::string> &seqs, 
    std::vector<std::pair <int, int> > &q_interval, 
    std::vector<std::pair <int, int> > &t_interval, 
    std::vector<int> &t_ID,
    std::vector<std::string> &q_seqs, 
    std::vector<std::string> &t_seqs
) {
  if(q_interval.size() != t_interval.size())  {
    cout << "Error: SequenceSearch::FetchAlignmentSeqs: inconsistent sizes of query and target intervals, abort." << endl;
    exit(1);
  }
  if(t_interval.size() != t_ID.size())  {
    cout << "Error: SequenceSearch::FetchAlignmentSeqs: inconsistent sizes of target intervals and target IDs, abort." << endl;
    exit(1);
  }
  int i;
  for(i = 0; i < q_interval.size(); ++ i) {
    q_seqs.push_back(
        query.substr(q_interval[i].first, q_interval[i].second - q_interval[i].first + 1)
    );
  }
  
  for(i = 0; i < t_interval.size(); ++ i) {
    t_seqs.push_back(
        seqs[t_ID[i]].substr(t_interval[i].first, t_interval[i].second - t_interval[i].first + 1)
    );
  }
  
  return;
}

void SequenceSearch::SelectTargetSeqs(
    const int mer_len, std::string &query, 
    std::vector<std::string> &candidate_seqs, 
    std::vector<std::string> &selected_seqs
) {
  if(query.length() < mer_len)  {
    cout << "Warning: SequenceSearch::SelectTargetSeqs: query sequence length too short!!!" << endl;
    return;
  }  
  int i, j;
  unordered_map<string, bool> q_mer;
  for(i = 0; i <= query.length() - mer_len; ++ i) {
    q_mer[query.substr(i, mer_len)] = true;
  }
  for(i = 0; i < candidate_seqs.size(); ++ i) {
    int num_hits = 0;
    for(j = 0; j <= candidate_seqs[i].length() - mer_len; ++ j) {
      if(q_mer.find(candidate_seqs[i].substr(j, mer_len)) != q_mer.end()) {
        ++ num_hits;
        if(num_hits >= 2) {
          selected_seqs.push_back(candidate_seqs[i]);
          continue;
        }
      }
    }
  }
  return;
}
