#include "read_cluster.h"

using namespace std;

ReadCluster::ReadCluster(
    const int& n_back, const std::string& file_name, const std::string& lcp_file, 
    const std::string& mcp_file, const std::string& gsa_file
)  {
  n_back_len_ = n_back;
  LoadSequence(file_name);
  LoadSuffixArray(lcp_file, mcp_file, gsa_file);
  reads_resolved = vector<bool> (suffix_array_->getSize(), false);
  return;
}

ReadCluster::~ReadCluster() {
  return;
}

void ReadCluster::LoadSequence(const std::string& file_name)  {
  vector<string> files_in;
	files_in.push_back(file_name);
	num_sequences_ = (unsigned int) seq::totalSequenceCount(files_in);
	sequence_ = new char *[num_sequences_];
  char **tag_foo = NULL;
  seq::loadSequences(files_in, tag_foo, sequence_, SEQONLY);
  delete tag_foo;
  return;
}

void ReadCluster::LoadSuffixArray(
    const std::string& lcp_file, const std::string& mcp_file, const std::string& gsa_file
) {
  suffix_array_ = new GSA();
  suffix_array_->load(
    lcp_file.c_str(), mcp_file.c_str(), gsa_file.c_str()
  );
  suffix_array_->setSequences(sequence_);
  suffix_array_->setReadCount(num_sequences_);
  return;
}

void ReadCluster::RecruitConnectedReads(void) {
  IndexSample encoder(n_back_len_, ALL20);
  unordered_map<int, int> read_partitioner;
  int part_id = 0, range_begin = 0, num_additional_paths = 0;
  bool set_begin = false;
  for(int i = 1; i < suffix_array_->getSize(); ++ i) {
    //if(i % 1000000 == 0)  {
    //  cout << i << " positions visited:  size: " << kmer_rank_.size()  << endl;
   // }
    if(set_begin && 
        ((int) suffix_array_->getLcp(i) < n_back_len_ || i == suffix_array_->getSize() - 1 || reads_resolved[i])
    )  {
      // defined a range, go through the range to identify unique kmers
      if(num_additional_paths == 0)  {
        string last_seq = suffix_array_->getSuffix(i - 1);
        string fw_seq = last_seq.substr(0, n_back_len_);
        string re_seq = last_seq.substr(last_seq.length() - n_back_len_, n_back_len_);
        KmerType fw_encode = encoder.encode_kmer(fw_seq);
        KmerType re_encode = encoder.encode_kmer(re_seq);
        int fw_rank, re_rank;
        auto it_fw = kmer_rank_.find(fw_encode);
        if(it_fw == kmer_rank_.end())  {
          pair<int, int> fw_bound = suffix_array_->searchWithLCPs((SfaChar*) fw_seq.c_str(), fw_seq.length());
          kmer_rank_[fw_encode] = fw_bound.first; 
          fw_rank = fw_bound.first;
        } else  {
          fw_rank = it_fw->second;
        }
        auto it_re = kmer_rank_.find(re_encode);
        if(it_re == kmer_rank_.end())  {
          pair<int, int> re_bound = suffix_array_->searchWithLCPs((SfaChar*) re_seq.c_str(), re_seq.length());
          kmer_rank_[re_encode] = re_bound.first;
          re_rank = re_bound.first;
        } else  {
          re_rank = it_re->second;
        }
        
        connected_kmers_[fw_rank][re_rank] = true;
        connected_kmers_[re_rank][fw_rank] = true;
      }
      set_begin = false;
    }
    if(set_begin && (int) suffix_array_->getLcp(i) < (int) suffix_array_->getSuffixLength(i - 1) && !reads_resolved[i])  {
      ++ num_additional_paths;
    }
    if(!set_begin && (int) suffix_array_->getSuffixLength(i) >= n_back_len_ && !reads_resolved[i])  {
      range_begin = i;
      num_additional_paths = 0;
      set_begin = true;
    }
  }
  return;
}

void ReadCluster::InterpretConnections(void)  {
  int index = 0;
  //cout << "size:  " << connected_kmers_.size() << endl;
  for(auto it = connected_kmers_.begin(); it != connected_kmers_.end(); ++ it)  {
    list<int> group_kmer;
    group_kmer.push_back(it->first);
    for(auto it_l = it->second.begin(); it_l != it->second.end(); ++ it_l)  {
      group_kmer.push_back(it_l->first);
    }
    // check if there is any overlap
    bool group_presented = false;
    int cluster_id = index;
    for(auto it_g = group_kmer.begin(); it_g != group_kmer.end(); ++ it_g) {
      auto it_recruited = kmer_cluster_.find(*it_g);
      if(it_recruited != kmer_cluster_.end()) {
        group_presented = true;
        cluster_id = it_recruited->second;
        break;
      }
    }
    // check if we need to initialize a new cluster
    if(!group_presented)  {
      ++ index;
    }
    // recruit the reads
    for(auto it_g = group_kmer.begin(); it_g != group_kmer.end(); ++ it_g)  {
      kmer_cluster_[*it_g] = cluster_id;    
    } 
  }
  //cout << "Num of clusters:  " << index << endl;
  /*
  while(!connected_kmers_.empty()) {
    cout << "Each round:  " << connected_kmers_.size() << endl; 
    queue<int> to_visit_kmers;
    unordered_map<int, bool> visited_kmers;
    to_visit_kmers.push(connected_kmers_.begin()->first);
    while(!to_visit_kmers.empty()) {
      int kmer = to_visit_kmers.front();
      to_visit_kmers.pop();
      if(visited_kmers.find(kmer) == visited_kmers.end())  {
        visited_kmers[kmer] = true;
        for(auto it = connected_kmers_[kmer].begin(); it != connected_kmers_[kmer].end(); ++ it) {
          to_visit_kmers.push(it->first);
        }
      }
    }
    for(auto it = visited_kmers.begin(); it != visited_kmers.end(); ++ it) {
      kmer_cluster_[it->first] = index;
      connected_kmers_.erase(it->first);
    }
    ++ index;
  }
  */
  return;
} 

void ReadCluster::GetClosureSequences(const std::list<std::string>& seed_kmers) {
  unordered_map<int, bool> presented_clusters;
  for(auto it = seed_kmers.begin(); it != seed_kmers.end(); ++ it) {
    pair<int, int> bound = suffix_array_->searchWithLCPs((SfaChar*) (*it).c_str(), (*it).length());
    int rank = bound.first;
    if(kmer_cluster_.find(rank) != kmer_cluster_.end())  {
      presented_clusters[kmer_cluster_[rank]] = true;
    } else  {
      recruited_kmer_[rank] = true;
    }
  }
  // define the recruited kmers
  for(auto it = kmer_cluster_.begin(); it != kmer_cluster_.end(); ++ it) {
    if(presented_clusters.find(it->second) != presented_clusters.end())  {
      recruited_kmer_[it->first] = true;
    }
  }
  // for each kmer define the reads
  unordered_map<int, bool> recruited_reads;
  for(auto it = recruited_kmer_.begin(); it != recruited_kmer_.end(); ++ it) {
    string kmer_sequence = suffix_array_->getSuffix(it->first);
    kmer_sequence = kmer_sequence.substr(0, n_back_len_);
    pair<int, int> bound = suffix_array_->searchWithLCPs((SfaChar*) kmer_sequence.c_str(), kmer_sequence.length());
    for(int i = bound.first; i <= bound.second; ++ i) {
      recruited_reads[suffix_array_->getId(i)] = true;
    }
  }
  // print the sequence
  for(auto it = recruited_reads.begin(); it != recruited_reads.end(); ++ it) {
    cout << ">read_" << it->first << endl;
    cout << suffix_array_->getSequence_explicit(it->first) << endl;
  }
  return;
}

std::string ReadCluster::GetReadSequence(const int& rank)  {
  return string(sequence_[rank]);
}
