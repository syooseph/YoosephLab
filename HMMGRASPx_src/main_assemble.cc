#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"
#include "reachable_reads.h"
#include "scoring_function.h"
#include "read_alignment.h"
#include "greedy_assembly.h"
#include "assemble_extend.h"
#include "contig_refinement.h"
#include "assemble_functor.h"

#include "sequence.h"
#include "timer.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

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
static int assembly_depth;
static double seed_score_scale = 0.7;
static int scoring_matrix = 0;
static int num_threads = 1;
static string verbose;

void PrintUsageBuild()  {
  cout << "Usage: graps-assemble query db contig_out" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
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
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value("WorkSpace"), "working directory for indexing file dump")
      ("seed_score_scale", boost::program_options::value<double>(&seed_score_scale), "selectivity for the seeds (default 0.6, higher means less seeds/faster runtime)")
      ("assembly_depth", boost::program_options::value<int>(&assembly_depth)->default_value(20), "maixmum depth of assembly path (each direction)")
      ("scoring_matrix", boost::program_options::value<int>(&scoring_matrix)->default_value(0), "scoring matrix\n 0: BLOSUM62, 1: BLOSUM80, 2: BLOSUM90,\n 3: BLOSUM50, 4: BLOSUM45, 5: PAM250,\n 6: PAM70, 7: PAM30")
      ("gap_open", boost::program_options::value<int>(&gap_open)->default_value(-11), "gap open penalty for sequence alignment")
      ("gap_extend", boost::program_options::value<int>(&gap_extend)->default_value(-1), "gap extension penalty for sequence alignment")
      ("band_size", boost::program_options::value<int>(&band_size)->default_value(20), "band size for sequence alignment")
      ("e_value", boost::program_options::value<double>(&e_value)->default_value(5.0), "e-value cutoff (Karlin-Altschul statistics)")
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
    cout << "Usage: graps-assemble query db out\n" << endl;
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
    cout << "Error: grasp-assemble: db_file does not exist." << endl;
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
  double start_time = mytime();
  double check_time;
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASPx::Assemble info: Begin of program execution." << endl;
  }
  // load sequences
  SequenceBuild* db_seq = new SequenceBuild(db_file);
  string db_stem = db_seq->GetFileStem(db_file);
  //db_seq.LoadSFA(workspace_dir, db_stem);
  string query_stem = db_seq->GetFileStem(query);
  // initiate scoreing function
  ScoringFunction<int>* score_scheme = new ScoringFunction<int>(
      PROTEIN, (MatrixName) scoring_matrix, gap_extend, gap_open
  );
  // select reads for seeding purpose
  string idx_map_file = workspace_dir + "/" + db_stem + ".rdm";
  string idx_seed_file = workspace_dir + "/" + db_stem + ".sxt"; 
  string idx_read_file = workspace_dir + "/" + db_stem + ".rxt";
  string idx_hsm_file = workspace_dir + "/" + db_stem + ".hsm";
  ReachableReads *candidate_selector = new ReachableReads(
      idx_map_file, idx_seed_file, idx_read_file, idx_hsm_file
  );
  //***************************
  if(is_verbose)  {
    check_time = mytime();
    cout << "GRASPx::Assemble info: Load sequences and indexing files done. ";
    printElapsed(start_time, check_time, "");
    cout << "============================================================" << endl;
  }
  //***************************
  
  //vector<char> aa_alpha = {
  //    'P', 'G', 'E', 'K', 'R', 'Q', 'D', 'S', 'N', 'T', 
  //    'H', 'C', 'I', 'V', 'W', 'Y', 'F', 'A', 'L', 'M'
  //};
  //ReducedAlphabet reduc_alph(GBMR10);  
  //candidate_selector.ComputeHighScoreMatch(aa_alpha, 0.3, 3, score_scheme, reduc_alph);
  
  //double hs_time = mytime();
  //printElapsed(start_time, hs_time, "compute high-score match");
  //start_time = mytime();
  
  // load query sequences
  vector<string> files_in;
  files_in.resize(1);
  files_in[0] = query;
  int num_queries = (unsigned int) seq::totalSequenceCount(files_in);
  char **query_tags = new char *[num_queries];
  char **query_seqs = new char *[num_queries];
  seq::loadSequences(files_in, query_tags, query_seqs, TAGSEQ);
  
  // schedule the assembly workers
  boost::threadpool::pool assembly_worker_pool(num_threads);
  vector<list<ContigType>* >* contigs_holder = new vector<list<ContigType>* >;
  contigs_holder->resize(num_queries); 
  LockType *mutex_out = new LockType;
  LockType *mutex_seq = new LockType;
  LockType *mutex_score = new LockType;
  LockType *mutex_link = new LockType;
  for(int id = 0; id < num_queries; ++ id) {
    (*contigs_holder)[id] = new list<ContigType>;
    string q_tag = query_tags[id];
    string q_seq = query_seqs[id];
    //AssembleFunctor assemble_worker(
    //    0, q_seq, db_seq, score_scheme, candidate_selector, 
    //    seed_score_scale, collect_depth, band_size, e_value,
    //    (*contigs_holder)[id]
    //);
    //assemble_worker.Run();
    boost::shared_ptr<AssembleFunctor> job(
      new AssembleFunctor(
          id, is_verbose, 
          mutex_out, mutex_seq, mutex_score, mutex_link, 
          q_tag, q_seq, 
          db_seq, score_scheme, candidate_selector, 
          seed_score_scale, assembly_depth, band_size, e_value,
          (*contigs_holder)[id] 
      )
    );
    boost::threadpool::schedule(assembly_worker_pool, boost::bind(&AssembleFunctor::Run, job));
  }
  assembly_worker_pool.wait();
  if(is_verbose)  {
    cout << "============================================================" << endl;
  }
  // write assembled contigs
  start_time = mytime();
  AssembleExtend guided_assembly;
  for(int id = 0; id < num_queries; ++ id) {
    int n = 0;
    for(auto itc = (*contigs_holder)[id]->begin(); itc != (*contigs_holder)[id]->end(); ++ itc) { 
      stringstream outs;
      outs << query_tags[id] << "||contig_" << n;
      string barcode = outs.str();  
      string header = guided_assembly.EncodeHeader(barcode, *itc);
      //string header = "foo";
      out_fh << ">" << header << endl << itc->sequence << endl;
      ++ n;
    }
  }
  out_fh.close();
  if(is_verbose)  {
    check_time = mytime();
    cout << "GRASPx::Assemble info: write assembled contigs done. ";
    printElapsed(start_time, check_time, "");
    start_time = mytime();
  }
  
  // write alignment
  ofstream aln_out_fh;
  if(vm.count("alignment_out"))  {
    ContigRefinement recalibrate_contigs;
    aln_out_fh.open(aln_output.c_str(), ios_base::out);
    
    if(aln_out_fh.good())  {
      for(int id = 0; id < num_queries; ++ id) {
        string q_tag = query_tags[id];
        string print_content;
        recalibrate_contigs.FormatForPrint(q_tag, *((*contigs_holder)[id]), print_content);
        aln_out_fh << print_content;
      }
      aln_out_fh.close();
    } else  {
      cout << "Warning: GRASP unable to write alignment file" << endl;
    }
    if(is_verbose)  {
      check_time = mytime();
      cout << "GRASPx::Assemble info: write alignments done. ";
      printElapsed(start_time, check_time, "");
    }
  }
  
  if(is_verbose)  {
    cout << "============================================================" << endl;
  }
  // collect memories
  for(int id = 0; id < num_queries; ++ id) {
    delete (*contigs_holder)[id];
  }
  delete contigs_holder;
  for(int id = 0; id < num_queries; ++ id) {
    delete [] query_tags[id];
    delete [] query_seqs[id];
  }
  delete [] query_tags;
  delete [] query_seqs;
  delete mutex_out;
  delete mutex_seq;
  delete mutex_score;
  delete mutex_link;
  delete db_seq;
  delete score_scheme;
  delete candidate_selector;
  return 0;
}
