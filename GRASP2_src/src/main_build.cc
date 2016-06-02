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

#include <iostream>
#include <string>
#include <list>
#include <ctime>
#include <cmath>

using namespace std;

static string workspace_dir;
static string db_file;
static int alph_id;
static int extd_len;
static int num_threads = 1;
static string verbose;
static int scoring_matrix = 0;

static int mer_len = 3;

void PrintUsage()  {
  cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl;
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
      ("db_file", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("index", boost::program_options::value<string>(&workspace_dir)->default_value("index"), "working directory for indexing file dump")
      ("extension_len", boost::program_options::value<int>(&extd_len)->default_value(10), "minimum overlap length for path extension")
      ("alphabet", boost::program_options::value<int>(&alph_id)->default_value(4), "reduced alphabet to be used for seeding\n  0:ALL20  1:DSSP5  2:DSSP10  3:GBMR4\n  4:GBMR10  5:HSDM5  6:SDM6  7:MURPHY5\n  8:MURPHY10  9:TD5  10:TD10")
      ("num_threads", boost::program_options::value<int>(&num_threads)->default_value(1), "maximum number of threads to be used")
      ("verbose", boost::program_options::value<string>(&verbose), "print intermediate information (default true)")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("db_file", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: grasp-build [peptide_db (FASTA)]" << endl << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  boost::filesystem::path abs_workspace = workspace_dir;
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-build: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::is_directory(workspace_dir))  {
    cout << workspace_dir << endl;
    cout << "Error: grasp-build: working space does not exist (please provide full path)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(extd_len < 6 || extd_len > 20)  {
    cout << "Error: grasp-build: extension length out of range (allowed range: 6-20)." << endl;
    exit(0);
  }
  if(alph_id < 0 || alph_id > 10)  {
    cout << "Error: grasp-build: reduced alphabet not supported (supported alphabets: 0-10)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  bool is_verbose = true;
  if(verbose == "False" || verbose == "false" || verbose == "No" || verbose == "no" || verbose == "0")  {
    is_verbose = false;
  }
  
  
  if(is_verbose)  {
    cout << "============================================================" << endl;
    cout << "GRASP2-Build: Begin of program execution." << endl;
  }
  
  BioAlphabet protein_alphabet(PROT);
  ReducedAlphabet reduced_alphabet((enum Alphabet) alph_id);
  ScoringProt scoring_function(static_cast<enum MatrixName>(scoring_matrix), -10, -1); 
  
  // Load in the peptide sequences to be searched against
  double start_time = MyTime();
  double check_time;
  Loader pepdb_loader;
  int num_seqs = pepdb_loader.CountFastaNumSeqs(db_file.c_str());
  char **header = new char* [num_seqs];
  char **seqs = new char* [num_seqs];
  num_seqs = pepdb_loader.LoadFasta(protein_alphabet, db_file.c_str(), header, seqs);
  string concat_seq; 
  Concatenator concat_obj(seqs, num_seqs, concat_seq);
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Load peptide database done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Construct Burrows-Wheeler Transformation
  start_time = MyTime();
  BWT bwt;
  bwt.ConstructNoPos(protein_alphabet, concat_seq.c_str());  
  concat_seq = string(concat_seq.rbegin(), concat_seq.rend());
  BWT rev_bwt;
  rev_bwt.ConstructNoPos(protein_alphabet, concat_seq.c_str());
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct Burrows-Wheeler transformation done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Compute pariwse overlap
  start_time = MyTime();
  StringGraph strG;
  std::vector<ExtType> extension;
  strG.MultiComputeExtension(
      num_threads, extd_len, 0, num_seqs, seqs, 
      bwt, rev_bwt, extension
  );
  strG.ImportExtension(extension);
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Construct string graph done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Post-processing of the string graph
  start_time = MyTime();
  strG.CheckGraph(); 
  while(strG.RemoveTipsBeforeCondense())  {;}
  strG.CondenseGraph(seqs);
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Post-process string graph done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Write the string graph to hard dist
  start_time = MyTime();
  string db_stem = GetFileStem(db_file);
  string idx_unitig_file = workspace_dir + "/" + db_stem + ".utg";  
  strG.WriteGraph(protein_alphabet, seqs, idx_unitig_file);
  strG.Purge();
  // Collect memory
  for(int idm = 0; idm < num_seqs; ++ idm) {
    delete [] header[idm]; delete [] seqs[idm]; 
  }
  delete [] header; delete [] seqs;
  //***************************
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Write string graph unitigs done.";
    PrintElapsed(start_time, check_time, "");
  }
  //***************************
  // Constructing k-mer mapping information
  //start_time = MyTime();
  //StringGraph strLoad;
  //vector<int> orphan_id;
  ////vector<string> orphan_seq;
  //strLoad.LoadGraph(idx_unitig_file, orphan_id, orphan_seq);
  //vector<string> unitigs;
  //strLoad.RecordEdgeSeqs(unitigs);
  //strLoad.Purge();
  //SequenceSearch seq_search; 
  //string idx_kmer_file = workspace_dir + "/" + db_stem + ".kmp"; 
  // include the single sequences into the indexing step
  //for(int i = 0; i < orphan_seq.size(); ++ i) {
  //  unitigs.push_back(orphan_seq[i]);
  //}
  //seq_search.IndexKmerPosition(protein_alphabet, mer_len, unitigs, idx_kmer_file);
  //***************************
  //if(is_verbose)  {
  //  check_time = MyTime();
  //  cout << "GRASP2-Build: Constructing and writing k-mer positions done. ";
  //  PrintElapsed(start_time, check_time, "");
    //cout << "GRASP2-Build: End of program execution." << endl;
    //cout << "============================================================" << endl;
  //}
  //***************************
  
  start_time = MyTime();
  SequenceSearch seq_search; 
  string idx_neighbor_file = workspace_dir + "/" + db_stem + ".knb";
  seq_search.IndexKmerNeighbor(
      mer_len, protein_alphabet, scoring_function, 
      11, idx_neighbor_file
  );  
  if(is_verbose)  {
    check_time = MyTime();
    cout << "GRASP2-Build: Constructing and writing k-mer neighbors done. ";
    PrintElapsed(start_time, check_time, "");
    cout << "GRASP2-Build: End of program execution." << endl;
    cout << "============================================================" << endl;
  }
  
  //***************************
  return 0;
}
