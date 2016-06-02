#include "../include/iterative_assembly.h"

using namespace std;

void PrintTime()  {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
  std::cout << "currentDateTime()=" << buf << std::endl;
  return;
}

bool IterativeAssembly::CheckConfig() {
  return true;
}

void IterativeAssembly::Assemble(
    std::vector<UnitigType> &unitigs, std::vector<bool> &read_mark
) {
  // initilize the read mark, assume all reads are not mapped
  read_mark.assign(num_seqs_, false);
  int mer_len = 6;
  KmerUnitcoder encoder(alphabet_, mer_len);
  BWTSearch bwt_searcher;
  FrequencyTable<KmerUnitType> unit_freq;
  ConstructFreqTable(mer_len, read_mark, unit_freq);
  // assemble iteratively
  int i;
  for(i = 0; i < kmer_config_.size(), i < coverage_config_.size(); ++ i) {
    // computing filtering frequency table
    UpdateFreqTable(mer_len, read_mark, unit_freq);
    cout << "Size of the frequency table: " << unit_freq.GetSize() << endl;
    //unit_freq.PrintAll();
    cout << "Finish constructing frequency table" << endl;
    // construct the graph
    DeBruijnGraph db_graph(encoder, kmer_config_[i]);
    db_graph.Construct(num_seqs_, seqs_, read_mark, unit_freq, coverage_config_[i]);
    db_graph.PrintGraphSizes();
    cout << "Finish constructing de bruijn graph" << endl;
    // applying a set of graph pruning
    db_graph.RefineGraph(coverage_config_[i]);
    db_graph.FormMinSpanningTree(2);    
    db_graph.BreakBadCovBranch();
    db_graph.RemoveTipsHead();
    db_graph.RemoveTipsTail();
    db_graph.CondenseGraph();
    db_graph.BridgeThreading(10);
    
    cout << "Finish pruning de bruijn graph" << endl;
    // idenfy the unitigs
    list<string> unitigs;
    db_graph.TraverseUnbranched(unitigs);
    cout << "Finish identifying unitigs" << endl;
    
#ifdef DEBUG
    int cid = 0;
    for(auto it = unitigs.begin(); it != unitigs.end(); ++ it) {
      cout << ">unitig_" << cid << endl << *it << endl;
      ++ cid;
    }
#endif
    
    PrintTime();
    cout << unitigs.size() << " unitigs are identified in this stage." << endl;
    
    // constructing BWT from the identified unitigs
    string concat_seq;
    Concatenator concat_obj(unitigs, concat_seq);  
    BWT bwt;
    bwt.Construct(alphabet_, concat_seq.c_str());  
    
    // TODO: add in-place reverse of concat_seq
    //concat_seq = string(concat_seq.rbegin(), concat_seq.rend());
    //concat_seq = concat_seq.substr(1, concat_seq.length() - 1) + concat_seq[0];
    concat_seq = string(concat_seq.rbegin(), concat_seq.rend());
    BWT rev_bwt;
    rev_bwt.Construct(alphabet_, concat_seq.c_str());
    cout << "Finish constructnig BWT" << endl;
    
    // align the reads to the unitigs
    int num_reads_taken = 0;
    AlignType pos;
    for(int i = 0; i < num_seqs_; ++ i) {
      if(read_mark[i]) continue;
      bwt_searcher.Search(bwt, rev_bwt, seqs_[i], pos);
      //cout << seqs_[i] << endl;
      // mark the read as taken
      if(strlen(seqs_[i]) - pos.q_pos >= 20) {
        read_mark[i] = true;
        ++ num_reads_taken;
      }
    }
    cout << "Finish mapping reads" << endl;
    
    PrintTime();
    cout << num_reads_taken << " reads are included in this stage." << endl;
    
    //for(int i = 0; i < num_seqs_; ++ i) {
    //  if(!read_mark[i])
    //    cout << ">read_" << i << endl << seqs_[i] << endl;
    //}
    //exit(0);
  }
  return;
}

void IterativeAssembly::ConstructFreqTable(
    const int mer_len, 
    std::vector<bool> &read_mark,
    FrequencyTable<KmerUnitType> &unit_freq
) {
  KmerUnitcoder unitcoder(alphabet_, mer_len);
  for(int i = 0; i < num_seqs_; ++ i) {
    if(read_mark[i]) continue;
    KmerUnitType c;
    for(int j = 0; j < strlen(seqs_[i]) - mer_len; ++ j) {
      if(j == 0) c = unitcoder.Encode(&seqs_[i][0]);
      else c = unitcoder.RightExt(c, seqs_[i][j + mer_len - 1]);
      unit_freq.Increase(c);
    }
  }
  return;
}

void IterativeAssembly::UpdateFreqTable(
    const int mer_len, 
    std::vector<bool> &read_mark,
    FrequencyTable<KmerUnitType> &unit_freq
) {
  KmerUnitcoder unitcoder(alphabet_, mer_len);
  for(int i = 0; i < num_seqs_; ++ i) {
    if(!read_mark[i]) continue;
    KmerUnitType c;
    for(int j = 0; j < strlen(seqs_[i]) - mer_len; ++ j) {
      if(j == 0) c = unitcoder.Encode(&seqs_[i][0]);
      else c = unitcoder.RightExt(c, seqs_[i][j + mer_len - 1]);
      unit_freq.Decrease(c);
    }
  }
  return;
}
