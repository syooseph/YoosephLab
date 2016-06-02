#include "../include/align_batch.h"
#include "../include/bwt.h"
#include "../include/bwt_search.h"
#include "../include/loader.h"
#include "../include/bio_alphabet.h"
#include "../include/kmer_unitcoder.h"
#include "../include/minimizer_sort.h"
#include "../include/string_graph.h"
#include "../include/sequence_search.h"
#include "../include/kmer_unitcoder.h"
#include "../include/scoring_prot.h"
#include "../include/reduced_alphabet.h"
#include "../include/kmer_filtering.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <omp.h>
#include <iostream>
#include <string>
#include <list>
#include <ctime>
#include <cmath>

using namespace std;

static string workspace_dir;
static string query;
static string db_file;
static string output;
static string aln_output;
static double e_value;
static int band_size;
static int gap_open;
static int gap_extend;
static double seed_score_scale = 0.9;
static int scoring_matrix = 0;
static int num_threads = 1;
static string verbose;

static int mer_len = 3;

void PrintUsage()  {
  cout << "Usage: grasp-assemble [query (FASTA)] [peptide_db (FASTA)] [contig_out]" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

double MyTime (void)
{
    int flag;
    clockid_t cid = CLOCK_REALTIME; // CLOCK_MONOTONE might be better
    timespec tp;
    double timing;
	
    flag = clock_gettime(cid, &tp);
    if (flag == 0) timing = tp.tv_sec + 1.0e-9*tp.tv_nsec;
    else           timing = -17.0;         // If timer failed, return non-valid time
	
    return(timing);
}

void PrintElapsed( double s, double e, const char *task )
{
	double elapsed = e - s ;
	printf ("[%02.0f:%02.0f:%05.2f]\t%s\n", 
			floor(elapsed/3600.0), 
			floor(fmod(elapsed,3600.0)/60.0), 
			fmod(elapsed,60.0),
			task);
	return;
}

string GetFileStem(const string& path)  {
  // going backward untill the '\/' character
  int i;
  for(i = path.length() - 1; i >= 0; -- i) {
    if(path[i] == '/')  break;
  }
  return path.substr(i + 1, path.length() - i - 1);
}

int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("query", boost::program_options::value<string>(&query), "query sequences (in FASTA format)")
      ("db", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("contig_out", boost::program_options::value<string>(&output), "assembled contigs output (in FASTA format)")
      ("alignment_out", boost::program_options::value<string>(&aln_output), "alignment output (between query and contigs, in BLAST-style format)")
      ("index", boost::program_options::value<string>(&workspace_dir)->default_value("index"), "working directory where the indexing files should be found")
      ("seed_score_scale", boost::program_options::value<double>(&seed_score_scale)->default_value(0.9), "selectivity for the seeds (between 0-1; higher value means less seeds/faster runtime)")
      ("scoring_matrix", boost::program_options::value<int>(&scoring_matrix)->default_value(0), "scoring matrix\n 0: BLOSUM62, 1: BLOSUM80, 2: BLOSUM90,\n 3: BLOSUM50, 4: BLOSUM45, 5: PAM250,\n 6: PAM70, 7: PAM30")
      ("gap_open", boost::program_options::value<int>(&gap_open)->default_value(-10), "gap open penalty for sequence alignment")
      ("gap_extend", boost::program_options::value<int>(&gap_extend)->default_value(-1), "gap extension penalty for sequence alignment")
      ("band_size", boost::program_options::value<int>(&band_size)->default_value(40), "band size for sequence alignment")
      ("e_value", boost::program_options::value<double>(&e_value)->default_value(10), "e-value cutoff (Karlin-Altschul statistics)")
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(1), "maximum number of threads to be used")
      ("verbose", boost::program_options::value<string>(&verbose), "print intermediate information (default true)")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("query", 1);
  pos_opt.add("db", 1);
  pos_opt.add("contig_out", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: grasp-assemble [query] [peptide_db] [contig_out]" << endl << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  if(!boost::filesystem::exists(query))  {
    cout << db_file << endl;
    cout << "Error: grasp-assemble: query does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-assemble: peptide database does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  ofstream out_fh;
  out_fh.open(output.c_str(), ios_base::out);
  if(!out_fh.good())  {
    cout << output << endl;
    cout << "Error: grasp-assemble: cannot write contig_out file." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::exists(workspace_dir))  {
    cout << "Error: grasp-assemble: working space does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  bool is_verbose = true;
  if(verbose == "False" || verbose == "false" || verbose == "No" || verbose == "no" || verbose == "0")  {
    is_verbose = false;
  }
  
  // Define basic alphabet and scoring function
  BioAlphabet *protein_alphabet = new BioAlphabet(PROT);
  ScoringProt *scoring_function = new ScoringProt(
      static_cast<enum MatrixName>(scoring_matrix), gap_open, gap_extend
  );
  AlignBatch *aligner = new AlignBatch;  

  //***************************
  double start_time = MyTime();
  double check_time;
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASP2-Assemble: Begin of program execution." << endl;
  }
  //***************************
  // Loading the query sequences
  Loader query_loader;
  int num_queries = query_loader.CountFastaNumSeqs(query.c_str());
  char **query_tags = new char *[num_queries];
  char **query_seqs = new char *[num_queries];
  query_loader.LoadFasta(*protein_alphabet, query.c_str(), query_tags, query_seqs);
  // Constructing the String Graph
  string db_stem = GetFileStem(db_file);
  StringGraph strG;
  vector<int> orphan_id;
  vector<string> orphan_seq;
  // ".UTG" stands for UniTiG
  string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  
  // ".RDM" stands for ReDuced Map
  string idx_kmer_file = workspace_dir + "/" + db_stem + ".kmp";  
  // ".KBF" stands for Kmer Bloom Filter  
  string idx_neighbor_file = workspace_dir + "/" + db_stem + ".knb";
  // Loading the unitig index file and extract the edge sequences
  strG.LoadGraph(idx_unitig_file, orphan_id, orphan_seq); 
  vector<string> *unitigs = new vector<string>; 
  int num_unitigs = strG.RecordEdgeSeqs(*unitigs);
  // include the single sequences into the indexing step
  for(int i = 0; i < orphan_seq.size(); ++ i) {
    unitigs->push_back(orphan_seq[i]);
  }
  vector<BoostSTREdge> graph_edge;
  strG.ComputeGraphEdgeMapping(graph_edge);
  // Compute the size of all unitigs
  long int db_size = 0;
  for(auto it = unitigs->begin(); it != unitigs->end(); ++ it) db_size += it->length();
  // Constructing the k-mer hashing table
  SequenceSearch seq_search;
  //seq_search.LoadKmerPosition(protein_alphabet, idx_kmer_file, mer_len, kmer_pos);
  vector<vector<int> > *kmer_neighbor = new vector<vector<int> >;
  seq_search.LoadKmerNeighbor(*protein_alphabet, idx_neighbor_file, mer_len, *kmer_neighbor);

  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Assemble: Load sequences and indexing files done. ";
    PrintElapsed(start_time, check_time, "");
    cout << "============================================================" << endl;
  }
  //*************************** 
  int **q_interval = new int* [num_queries];
  int **t_interval = new int* [num_queries];
  int **mx_score = new int* [num_queries];
  int **t_ID = new int* [num_queries]; 
  
  int q_id;
  #pragma omp parallel num_threads(num_threads) private(q_id) shared(protein_alphabet, scoring_function, unitigs, kmer_neighbor, q_interval, t_interval, mx_score, t_ID)
  {
    #pragma omp for
    for(q_id = 0; q_id < num_queries; ++ q_id) {
      double start_time_p, check_time_p;
      #pragma omp critical
      {
        start_time_p = MyTime();
      }
      string current_query = query_seqs[q_id];
      // Pre-processing and identify candidate unitigs for alignment filtering
      int screen_cutoff = scoring_function->ComputeCutoff(
          current_query.length(), db_size, e_value * 1000
      );
      int final_cutoff = scoring_function->ComputeCutoff(
          current_query.length(), db_size, e_value
      );
      //********************************
      
      q_interval[q_id] = new int [2 * unitigs->size()];
      t_interval[q_id] = new int [2 * unitigs->size()];
      mx_score[q_id] = new int [unitigs->size()];
      t_ID[q_id] = new int [unitigs->size()];
      
      int num_filtered = seq_search.SelectID(
          *protein_alphabet, mer_len, screen_cutoff, current_query, 
          *unitigs, *scoring_function, *kmer_neighbor, 
          q_interval[q_id], t_interval[q_id], mx_score[q_id], t_ID[q_id]
      );
      
      #pragma omp critical
      {   
        
        if(is_verbose)  {
          check_time_p = MyTime();
          cout << "GRASP2-Assemble::" << query_tags[q_id] << ": Unitig filtering done.";
          PrintElapsed(start_time_p, check_time_p, "");
          start_time_p = MyTime();
        }
      }
      //********************************

      // Constructing paths based on alignment scores of the unitigs
      vector<string> high_score_seqs;
      strG.GetHighScoringPaths(
          screen_cutoff, current_query.length(), *
          unitigs, graph_edge, 
          num_filtered, t_ID[q_id], mx_score[q_id], q_interval[q_id], 
          0.95, high_score_seqs
      );
      #pragma omp critical
      { 
        delete [] q_interval[q_id];
        delete [] t_interval[q_id];
        delete [] mx_score[q_id];
        delete [] t_ID[q_id];
      }
      #pragma omp critical
      {  
        if(is_verbose)  {
          check_time_p = MyTime();
          cout << "GRASP2-Assemble::" << query_tags[q_id] << ": Path reconstruction done.";
          PrintElapsed(start_time_p, check_time_p, "");
          start_time_p = MyTime();
        }
      }
      //***************************
      
      // Re-aligning the paths with the query     
      vector<int> score(high_score_seqs.size(), 0);
      for(int i = 0; i < high_score_seqs.size(); ++ i) {
        int b = band_size + 
            2 * abs((int) current_query.length() - 
                    (int) high_score_seqs[i].length()
            );
        score[i] = aligner->AlignLocalPairwise(
            current_query, high_score_seqs[i], 
            b, *scoring_function
        );
      }
      
      vector<string> out_seqs;
      for(int i = 0; i < score.size(); ++ i) {
        if(score[i] >= final_cutoff) 
          out_seqs.push_back(high_score_seqs[i]);
      }
      
      vector<FullAlnType> aln_info;
      aligner->MultiAlignLocalFull(
          1, current_query, out_seqs.size(), out_seqs, 
          *scoring_function, aln_info
      );

      //***************************
      #pragma omp critical
      {
        if(is_verbose)  {
          check_time_p = MyTime();
          cout << "GRASP2-Assemble::" << query_tags[q_id] << ": Re-alignment done.";
          PrintElapsed(start_time_p, check_time_p, "");
          start_time_p = MyTime();
        }
      }
      //***************************
      // Write assembled contigs
      #pragma omp critical
      {        
        for(int idx = 0; idx < aln_info.size(); ++ idx) {
          if(aln_info[idx].score < final_cutoff) continue;
          stringstream outs;
          outs << query_tags[q_id] << "||contig_" << idx;
          string barcode = outs.str();  
          int begin = aln_info[idx].t_interval.first;
          int end = aln_info[idx].t_interval.second;
          out_fh << ">" << barcode;
          out_fh << "||" << scoring_function->ComputeEValue(
              current_query.length(), db_size, aln_info[idx].score
          );
          out_fh << "||" << scoring_function->ComputeBitScore(aln_info[idx].score);
          out_fh << "||" << aln_info[idx].score << endl;
          out_fh << aln_info[idx].sequences.second.substr(begin, end - begin + 1) << endl;
        }
      }
      //***************************
      #pragma omp critical
      {
        if(is_verbose)  {
          check_time_p = MyTime();
          cout << "GRASP2-Assemble::" << query_tags[q_id] << ": Writing assembled contigs done.";
          PrintElapsed(start_time_p, check_time_p, "");
        }
      }
      //***************************
      
    }
    #pragma omp taskwait
  }
  out_fh.close();
  delete [] q_interval;
  delete [] t_interval;
  delete [] mx_score;
  delete [] t_ID;
  
  if(is_verbose)  {
    cout << "GRASP2-Assemble: End of program execution." << endl;
    cout << "============================================================" << endl;
  }
  // Collect memory
  for(int idm = 0; idm < num_queries; ++ idm) {
    delete [] query_tags[idm]; delete [] query_seqs[idm]; 
  }
  delete [] query_tags; delete [] query_seqs;
  return 0;
}
