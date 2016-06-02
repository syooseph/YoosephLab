#include "kmer_distance.h"

using namespace std;

KmerDistance::KmerDistance(void)  {
  return;
}

KmerDistance::KmerDistance(unsigned int in_kmer_size, enum Alphabet in_alphabet)
  : kmer_size_(in_kmer_size),
    alphabet_(in_alphabet)  
{
  return;
}

KmerDistance::~KmerDistance(void) {
  return;
} 

void KmerDistance::KmerAdjacency(string file_in)  {
  // open the file
  ifstream in_fh(file_in.c_str(), ios_base::in);
  assert(in_fh.good());
  // read the file
  int counter = 0;
  string line;
  string appended_seq = "";
  while(getline(in_fh, line))  {
    if(line[0] == '>')  {
      // update the kmer adjacency list with the sequence
      if(appended_seq != "")  {
        UpdateAdjacency(appended_seq);
      }
      ++ counter;
      // clear the concaternated sequence
      appended_seq = "";
    }  else  {
      appended_seq += line;
    }
  }
  cerr << "Size of kmer vertices: " << kmer_adj_list_.size() << endl;
  return;
}

void KmerDistance::UpdateAdjacency(std::string seq) {
  int i, j;
  IndexSample encoder(kmer_size_, alphabet_);
  vector<string> encoded_kmers(seq.length() - kmer_size_ + 1);
  for(i = 0; i < encoded_kmers.size(); ++ i) {
    encoded_kmers[i] = seq.substr(i, kmer_size_);
  }
  // define the boundary kmers
  boundary_kmers_[encoded_kmers[0]] = boundary_kmers_[encoded_kmers.back()] = true;
  //for(i = 0; i < encoded_kmers.size() - 1; ++ i) {
  //  for(j = i + 1; j < encoded_kmers.size(); ++ j) {
  //    if(boundary_kmers_.find(encoded_kmers[i]) != boundary_kmers_.end() && 
  //        boundary_kmers_.find(encoded_kmers[j]) != boundary_kmers_.end())  {
        // assuming we only check the boundaries
        i = 0;
        j = encoded_kmers.size() - 1;
        kmer_adj_list_[encoded_kmers[i]][encoded_kmers[j]] = true;
        kmer_adj_list_[encoded_kmers[j]][encoded_kmers[i]] = true;
  //    }
  //  }
  //}
  //int total_num_edges = 0;
  //for(auto it = kmer_adj_list_.begin(); it != kmer_adj_list_.end(); ++ it)  {
  //  total_num_edges += (*it).second.size();
  //}
  //cerr << "Total number of edges:	" << (total_num_edges / 2)  << endl; 
  return;
}

void KmerDistance::GetKmerClosure(
    const double& conv_rate, 
    const std::unordered_map<string, bool>& seed_kmers, 
    std::unordered_map<string, bool>& kmer_closure
)  {
  int total_num_edges = 0;
  for(auto it = kmer_adj_list_.begin(); it != kmer_adj_list_.end(); ++ it)  {
    total_num_edges += (*it).second.size();
  }
  cerr << "Total number of edges:       " << (total_num_edges / 2)  << endl;
  unordered_map<string, bool> checked_kmer;
  //cerr << "Hello world" << endl;
  // define the initial candidate kmers
  for(auto it = seed_kmers.begin(); it != seed_kmers.end(); ++ it) {
    KmerIterType q_elem;
    q_elem.kmer = it->first;
    q_elem.iter = 0;
    candidate_kmer_.push(q_elem);
  }
  int prev_iter = 0;
  int prev_closure_size = 0; 
  //cerr << "candidate size:  " << candidate_kmer_.size() << endl;
  while(!candidate_kmer_.empty()) {
    KmerIterType q_elem = candidate_kmer_.front();
    candidate_kmer_.pop();
    //cerr << "kmers in iteration:  " << q_elem.iter << endl;
    if(q_elem.iter != prev_iter)  {
      // if the closure converges as defined, then return
      cerr << "kmer-closure size:	" << kmer_closure.size() << endl;
      if((double) prev_closure_size / (double) kmer_closure.size() > conv_rate)  {
        //cerr << "returned for convergence" << endl;
        return;
      }
      prev_closure_size = kmer_closure.size();
    }
    prev_iter = q_elem.iter;
    if(checked_kmer.find(q_elem.kmer) != checked_kmer.end())  {
      //cerr << "continued for presences" << endl;
      continue;
    } 
    // if the kmer is not present yet
    for(auto it = kmer_adj_list_[q_elem.kmer].begin(); it != kmer_adj_list_[q_elem.kmer].end(); ++ it) {
      if(it->second)  {
        kmer_closure[it->first] = true;
        //cerr << (unsigned int) it->first << endl;
        if(checked_kmer.find(it->first) == checked_kmer.end())  {
          //cerr << "Entered here" << endl;
          KmerIterType n_elem;
          n_elem.kmer = it->first;
          n_elem.iter = q_elem.iter + 1;
          candidate_kmer_.push(n_elem);
        }
      }
    }
    checked_kmer[q_elem.kmer] = true;
  }   
  return;
}


