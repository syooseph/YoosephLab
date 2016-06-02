#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"
#include "reachable_reads.h"
#include "assemble_extend_p.h"
#include "contig_refinement_p.h"
#include "assemble_functor_p.h"

#include "sequence.h"
#include "timer.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <list>
#include <map>

using namespace std;

static string workspace_dir;
static string query;
static string db_file;
static string output;
static string aln_output;
static double p_value = 0.05;
static int band_size;
static int assembly_depth;
static double seed_score_scale = 0.8;
static int num_threads = 1;
static string verbose;
static int progressive_extend = 1;
static string dup_reads;

void PrintUsageBuild()  {
  cout << "Usage: graspxp-assemble query db contig_out" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

void ProcessMultiHMM(std::string &m_hmm_file, std::list<HMMProfile> &queries)  {
  ifstream HMMFILE(m_hmm_file);
  if(!HMMFILE.is_open())  {
    cout << "Main::Error: Cannot open HMM profile file." << endl;
    exit(1);
  }
  // reads in the HMM info
  std::list<std::string> hmm_content;
  std::string line;
  while(getline(HMMFILE, line)) {
    hmm_content.push_back(line);
    if(line == "//")  {
      HMMProfile q(hmm_content);
      queries.push_back(q);
      hmm_content.clear();
    }
  }
  HMMFILE.close();
  return;
}

bool _cmp_model(HMMProfile &a, HMMProfile &b) {
  if(a.GetProfileLen() >= b.GetProfileLen())  {
    return true;
  }
  return false;
}

int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("query", boost::program_options::value<string>(&query), "query profile (in HMMER3 hmm profile model format)")
      ("db", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("contig_out", boost::program_options::value<string>(&output), "assembled contigs output (in FASTA format)")
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value("WorkSpace"), "working directory for indexing file dump")
      ("seed_score_scale", boost::program_options::value<double>(&seed_score_scale), "selectivity for the seeds (default 0.7, higher means less seeds/faster runtime)")
      ("assembly_depth", boost::program_options::value<int>(&assembly_depth)->default_value(5), "maixmum depth of assembly path (each direction)")
      ("band_size", boost::program_options::value<int>(&band_size)->default_value(20), "band size for sequence alignment")
      ("p_value", boost::program_options::value<double>(&p_value), "p-value cutoff, default 0.05")
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(1), "maximum number of threads to be used")
      ("dup_reads", boost::program_options::value<string>(&dup_reads), "allowing seeds to be used multiple times for higher sensitivity (default false, --dup_reads=true/false)")
      ("progressive_extend", boost::program_options::value<int>(&progressive_extend), "number of steps to extend un-branched contig assembly (default: 1)")
      ("verbose", boost::program_options::value<string>(&verbose), "print intermediate information (default true, --verbose=true/false)")
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
    cout << "Usage: graspxp-assemble query db out\n" << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  if(!boost::filesystem::exists(query))  {
    cout << db_file << endl;
    cout << "Error: graspxp-assemble: query profile does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: graspxp-assemble: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  std::ofstream out_fh;
  out_fh.open(output.c_str(), std::ios_base::out);
  if(!out_fh.good())  {
    std::cout << output << std::endl;
    std::cout << "Error: graspxp-assemble: cannot write contig_out file." << std::endl;
    std::cout << "Please use \'--help\' for more details." << std::endl;
    exit(0);
  }
  out_fh.close();
  if(!boost::filesystem::exists(workspace_dir))  {
    cout << "Error: graspxp-assemble: working space does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  bool is_verbose = true;
  if(verbose == "False" || verbose == "false" || verbose == "No" || verbose == "no" || verbose == "0")  {
    is_verbose = false;
  }
  bool is_dup_reads= false;
  if(dup_reads == "True" || dup_reads == "true" || dup_reads == "Yes" || dup_reads == "yes" || dup_reads == "1")  {
    is_dup_reads = true;
  }
  double start_time = mytime();
  double check_time;
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASPxp::Assemble info: Begin of program execution." << endl;
  }
  // load sequences
  SequenceBuild* db_seq = new SequenceBuild(db_file);
  string db_stem = db_seq->GetFileStem(db_file);
  string query_stem = db_seq->GetFileStem(query);
  // initiate scoreing function
  // select reads for seeding purpose
  string idx_map_file = workspace_dir + "/" + db_stem + ".rdm";
  string idx_seed_file = workspace_dir + "/" + db_stem + ".sxt"; 
  string idx_read_file = workspace_dir + "/" + db_stem + ".rxt";
  string idx_hsm_file = workspace_dir + "/" + db_stem + ".hsm";
  ReachableReads *candidate_selector = new ReachableReads(
      *db_seq,
      idx_map_file, idx_seed_file, idx_read_file, idx_hsm_file
  );
  candidate_selector->ReadExtToVectors();
  //***************************
  if(is_verbose)  {
    check_time = mytime();
    cout << "GRASPxp::Assemble info: Load sequences and indexing files done. ";
    printElapsed(start_time, check_time, "");
    cout << "============================================================" << endl;
  }
  //***************************
  start_time = mytime();
  // load query profile
  list<HMMProfile> models;
  ProcessMultiHMM(query, models);
  models.sort(_cmp_model);
  
  boost::threadpool::pool assembly_worker_pool(num_threads);
  LockType *mutex_out = new LockType;
  LockType *mutex_seq = new LockType;
  LockType *mutex_score = new LockType;
  LockType *mutex_link = new LockType;
  int id = 0;
  for(auto it_m = models.begin(); it_m != models.end(); ++ it_m) {
    boost::shared_ptr<AssembleFunctorP> job(
      new AssembleFunctorP(
          id, is_verbose, 
          mutex_out, mutex_seq, mutex_score, mutex_link, 
          *it_m, db_seq, candidate_selector, 
          seed_score_scale, progressive_extend, is_dup_reads, 
          assembly_depth, band_size, p_value,
          output 
      )
    );
    boost::threadpool::schedule(assembly_worker_pool, boost::bind(&AssembleFunctorP::Run, job));
    ++ id;
  }
  assembly_worker_pool.wait();
  if(is_verbose)  {
    cout << "============================================================" << endl;
  }


  delete mutex_out;
  delete mutex_seq;
  delete mutex_score;
  delete mutex_link;
  delete db_seq;
  delete candidate_selector;
    
  return 0;
}
