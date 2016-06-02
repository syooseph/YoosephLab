#ifndef _SEQUENCE_SEARCH_H_
#define _SEQUENCE_SEARCH_H_

#ifndef MAX
#define MAX 999999
#endif

// the class is used to perform single-sequence-based search against the read set

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <tuple>

#include "align_batch.h"
#include "scoring_prot.h"
#include "bio_alphabet.h"
#include "kmer_unitcoder.h"
#include "kmer_filtering.h"
#include "reduced_alphabet.h"

struct AlignIntervalType  {
  // q1, q2 are locations in the query sequence
  // t1, t2 are locations in the target sequence
  int q1, q2, t1, t2;
  // score is the goodness of the alignment
  int score;
};



class SequenceSearch  {
 public:
  explicit SequenceSearch(void) {}
  ~SequenceSearch() {}
  
  // extract the k-mers that present in all "seqs"
  // the output "reduced_map" is a hash table, the key is the string in the reduced alphabet
  // the value is a list of k-mers that can be reduced into the corresponding key value
  // the output "seed_mer" is a hash table, the key is different k-mers
  // the value is a list of edge sequences that contain the corresponding k-mer
  void ExtractSeedMers(
      const int mer_len, ReducedAlphabet &alpha, std::vector<std::string> &seqs,
      std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
      std::unordered_map<std::string, std::vector<int> > &seed_mer
  );
  
  void ComputeReducedMap(
      const int mer_len, BioAlphabet &alpha,
      ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
      std::unordered_map<std::string, std::set<KmerUnitType> > &reduced_map
  );
  
  void WriteReducedMap(
      const std::string &file, const int mer_len, BioAlphabet &alpha,
      std::unordered_map<std::string, std::set<KmerUnitType> > &reduced_map
  );
  
  void IndexReducedMap(
      const int mer_len, BioAlphabet &alpha,
      ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
      const std::string &out_file
  );
  
  void LoadReducedMap(
      const std::string &in_file, 
      std::unordered_map<std::string, std::vector<std::string> > &reduced_map
  );
  
  void ComputeSeedMap(
      const int mer_len, ReducedAlphabet &alpha, std::vector<std::string> &seqs,
      std::unordered_map<std::string, std::vector<int> > &seed_mer
  );
  
  
  void WriteSeedMap(
      const std::string &file, 
      std::unordered_map<std::string, std::vector<int> > &seed_mer
  );
  
  void ComputeKmerFilter(
      const int mer_len, BioAlphabet &alpha, std::vector<std::string> &seqs,
      std::vector<std::string> &all_filters  
  );
  
  void WriteKmerFilter(const std::string &file, std::vector<std::string> &all_filters);
  //************************  
  void IndexKmerFilter(
      const int mer_len, BioAlphabet &alpha,
      ReducedAlphabet &re_alpha, std::vector<std::string> &seqs,
      const std::string &out_file
  );
  void LoadKmerFilter(const std::string &file, std::vector<KmerFiltering> &bm_filters);  
  //************************  
  void IndexKmerPosition(
      BioAlphabet &alphabet, const int mer_len, 
      std::vector<std::string> &seqs, const std::string &out_file
  );
  void LoadKmerPosition(
      BioAlphabet &alphabet, const std::string &in_file,
      const int mer_len, std::vector<std::vector<int> >& kmer_pos
  ); 
  //************************  
  void IndexKmerNeighbor(
      const int mer_len, BioAlphabet &alpha, 
      ScoringProt &fs, const int neighbor_score,
      const std::string &out_file
  );
  void LoadKmerNeighbor(
      BioAlphabet &alphabet, const std::string &in_file, const int mer_len,
      std::vector<std::vector<int> >& kmer_neighbor
  );
  //************************  
  
  void IndexKmerTarget(
      const int mer_len, BioAlphabet &alpha,
      std::vector<std::string> &seqs, 
      ScoringProt &fs, const int neighbor_score,
      const std::string &out_file
  );
  
  // write the seed-mers to hard-disk for future loading
  void DumpSeedMers(
      const std::string &file,
      std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
      std::unordered_map<std::string, std::vector<int> > &seed_mer
  );
  
  // load the seed-mers from hard-dist
  void LoadSeedMers(
      const std::string &file,
      std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
      std::unordered_map<std::string, std::vector<int> > &seed_mer
  );
  
  // given the query, compute the set of sequences that share a high-scoring kmer
  // with the query; "fs" is the scoring function, "scale" is the scale to define
  // high-scoring (any k-mer match with a score > scale * fs_diagonal * mer_len)
  // is considered as high-scoring; resulting high-score read IDs are stored in selected_IDs
  // the key of "selected_ID" is the edge sequence ID; the value is a list of seed pairs;
  // where the first value is the kmer and the second value is the location in the query sequence
  void SelectID(
      BioAlphabet &alphabet, const int mer_len, 
      const int screen_score, std::string &query, 
      std::vector<std::string> &seqs, ScoringProt &fs,
      std::vector<std::vector<int> >& kmer_neighbor,
      std::vector<std::pair <int, int> > &q_interval, 
      std::vector<std::pair <int, int> > &t_interval,
      std::vector<int> &mx_score, std::vector<int> &t_ID
  );
  
  int SelectID(
      BioAlphabet &alphabet, const int mer_len, 
      const int screen_score, std::string &query, 
      std::vector<std::string> &seqs, ScoringProt &fs,
      std::vector<std::vector<int> >& kmer_neighbor,
      int *q_interval, int *t_interval,
      int *mx_score, int *t_ID
  );
  
  void SelectID(
      const int mer_len, BioAlphabet &alpha, ReducedAlphabet &re_alpha, 
      std::string &query, std::vector<std::string> &seqs,
      std::unordered_map<std::string, std::vector<std::string> > &reduced_map,
      std::vector<KmerFiltering> &bm_filters, double scale, ScoringProt &fs, 
      std::vector<std::vector<int> > &seed_positions
  );
  
  // given the "query" and "target" sequences, and their kmer (with length "mer_len")
  // matching locations, find the interval in the query ("q_interval") 
  // and the interval in the target ("t_interval") where sequence should be extracted
  // returns the number of matchings taken
  int FindBestMatching(
      const int mer_len, const int screen_score,
      const std::string& query, const std::string &target, 
      const std::vector<AlignIntervalType> &matching, ScoringProt &fs, 
      std::vector<std::pair<int, int> > &q_interval, 
      std::vector<std::pair<int, int> > &t_interval,
      std::vector<int> &mx_pair_score
  );
  
  void DefineChainingPreSet(
      const int band_size,
      const std::vector<AlignIntervalType> &matching, 
      std::vector<std::vector<int> > &pre_set
  );
  
  // given the partial seeding information "selected_ID", define a set of non-overlapping regions
  // to perform the actual alignment. "mer_len" is the length of the kmer, "query" is the query sequence
  // "seqs" is the set of target sequences (edges in the graph), "selected_ID" is a set of 
  // target sequences that share common seeds with the query, "q_interval" is a set of query
  // intervals to be aligned, "t_interval" is a set of target intervals to be aligned,
  // "t_ID" is a set of IDs that maps the corresponding alignment to indexes in "seqs"
  void SelectAlignmentRegions(
      const int mer_len, const int screen_score, ScoringProt &fs, 
      std::string &query, std::vector<std::string> &seqs, 
      std::vector<std::vector<int> > &seed_positions,
      std::vector<std::pair <int, int> > &q_interval, 
      std::vector<std::pair <int, int> > &t_interval,
      std::vector<int> &mx_score,
      std::vector<int> &t_ID    
  );
  
  // copy the sequences to be aligned into "q_seqs" and "t_seqs"
  // "query" is the query sequence, "seqs" is the target sequences
  // "q_interval" and "t_intervals" are intervals in the query and target, respectively
  // "t_ID" is the collection of target IDs
  void FetchAlignmentSeqs(
      std::string &query, std::vector<std::string> &seqs, 
      std::vector<std::pair <int, int> > &q_interval, 
      std::vector<std::pair <int, int> > &t_interval, 
      std::vector<int> &t_ID,
      std::vector<std::string> &q_seqs, 
      std::vector<std::string> &t_seqs
  );
  
  void SelectTargetSeqs(
      const int mer_len, std::string &query, 
      std::vector<std::string> &candidate_seqs, 
      std::vector<std::string> &selected_seqs
  );
  
  void ExtendSeedToMatching(
      const int mer_len, const std::string& query, const std::string &target, 
      const std::vector<int> &all_seeds, ScoringProt &fs,
      std::vector<AlignIntervalType> &matching 
  );
  
  void ExtendSeedToMatchingTwoHit(
      const int mer_len, const std::string& query, const std::string &target, 
      const std::vector<int> &all_seeds, ScoringProt &fs,
      std::vector<AlignIntervalType> &matching 
  );
  
 private: 
};


#endif
