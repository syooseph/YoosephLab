#include "sequence_build.h"
#include "reduced_alphabet.h"
#include "database_index.h"

#include "timer.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

using namespace std;

static string workspace_dir;
static string db_file;
static int seed_len;
static int extd_len;
static int alph_id;
static int min_seed_coverage;
static int min_ext_coverage;

void PrintUsageBuild()  {
  cout << "Usage: graps-build db_file" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("db_file", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value("WorkSpace"), "working directory for indexing file dump")
      ("seed_coverage", boost::program_options::value<int>(&min_seed_coverage)->default_value(3), "minimum coverage required for seed sequence")
      ("extension_coverage", boost::program_options::value<int>(&min_ext_coverage)->default_value(1), "minimum coverage required for extending assembly")
      ("seed_len", boost::program_options::value<int>(&seed_len)->default_value(6), "length of the seeds")
      ("extension_len", boost::program_options::value<int>(&extd_len)->default_value(10), "minimum overlap length for path extension")
      ("alphabet", boost::program_options::value<int>(&alph_id)->default_value(4), "reduced alphabet to be used for seeding\n  0:ALL20  1:DSSP5  2:DSSP10  3:GBMR4\n  4:GBMR10  5:HSDM5  6:SDM6  7:MURPHY5\n  8:MURPHY10  9:TD5  10:TD10")
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
    cout << "Usage: graps-build db_file(expecting FASTA-format)\n" << endl;
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
  if(seed_len < 3 || seed_len > 10)  {
    cout << "Error: grasp-build: seed length out of range (allowed range: 3-10)." << endl;
    exit(0);
  }
  if(extd_len < 6 || seed_len > 20)  {
    cout << "Error: grasp-build: extension length out of range (allowed range: 6-20)." << endl;
    exit(0);
  }
  if(alph_id < 0 || alph_id > 10)  {
    cout << "Error: grasp-build: reduced alphabet not supported (supported alphabets: 0-10)." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  cout << "============================================================" << endl;
  double start_time = mytime();
  // load in sequence from file for forward sequences
  SequenceBuild db_seq(db_file);
  string db_stem = db_seq.GetFileStem(db_file);
  db_seq.BuildSFADefault();
  db_seq.BuildKeyArrayDefault();
  double fw_sfa_build_time = mytime();
  printElapsed(start_time, fw_sfa_build_time, "GRASPx::Build info: finished building forward suffix array");
  start_time = mytime();
  // dump the forward version suffix array
  db_seq.DumpSFA(workspace_dir, db_stem);
  double fw_sfa_dump_time = mytime();
  printElapsed(start_time, fw_sfa_dump_time, "GRASPx::Build info: finished writing forward suffix array");
  start_time = mytime();
  // prepare DatabaseIndex object
  DatabaseIndex db_dump(alph_id, seed_len, extd_len);
  db_dump.BuildSeedmerMap(db_seq);
  unordered_map<std::string, std::list<std::string> > reduc_alph_map;
  db_dump.CreateReducedMap(reduc_alph_map);
  string out_file = workspace_dir + "/" + db_stem + ".rdm";
  db_dump.DumpReducedMap(reduc_alph_map, out_file);
  // create forward seed extensions
  unordered_map<std::string, std::list<PositionType> > fw_seed_ext;
  db_dump.GetSeedExt(db_seq, min_seed_coverage, fw_seed_ext);
  // create forward read extensions
  std::unordered_map<RIDType, std::list<OverlapType> > fw_read_ext;
  db_dump.CreateReadExtWorker(min_ext_coverage, db_seq, fw_read_ext);
  double fw_ext_time = mytime();
  printElapsed(start_time, fw_ext_time, "GRASPx::Build info: finished resolving forward extensions");
  start_time = mytime();
  // finish using forward suffix array for now, destruct it to save memory
  db_seq.DestructSFA();
  // build the reverse suffix array first, as we can discarded immediately
  // after resolving extensions
  SequenceBuild db_seq_rev(db_seq);
  db_seq_rev.InPlaceReverse();
  // build reverse suffix array 
  db_seq_rev.BuildSFADefault();
  db_seq_rev.BuildKeyArrayDefault();
  double re_sfa_build_time = mytime();
  printElapsed(start_time, re_sfa_build_time, "GRASPx::Build info: finished building reverse suffix array");
  start_time = mytime();
  // create reverse seed_extensions
  unordered_map<std::string, std::list<PositionType> > re_seed_ext;
  db_dump.GetSeedExtRev(db_seq_rev, min_seed_coverage, re_seed_ext);
  // create reverse read extensions
  std::unordered_map<RIDType, std::list<OverlapType> > re_read_ext;
  db_dump.CreateReadExtWorker(min_ext_coverage, db_seq_rev, re_read_ext);
  db_seq_rev.DestructSFA();
  db_seq_rev.DestructSequences();
  double re_ext_time = mytime();
  printElapsed(start_time, re_ext_time, "GRASPx::Build info: finished resolving reverse extensions");
  start_time = mytime();
  // merge seed extension pairs and dump read extension
  // we need to load forward SFA to merge the seed paris
  db_seq.LoadSFA(workspace_dir, db_stem);
  unordered_map<string, list<ReadPairType> > seed_pair_ext;
  db_dump.MatchSeedPair(db_seq, fw_seed_ext, re_seed_ext, seed_pair_ext);
  double pair_time = mytime();
  printElapsed(start_time, pair_time, "GRASPx::Build info: finished bridging seed pairs");
  start_time = mytime();
  //db_dump.CreateSeedExt(min_seed_coverage, db_seq, db_seq_rev, seed_pair_ext);
  string out_file_seed_ext = workspace_dir + "/" + db_stem + ".sxt";
  db_dump.DumpSeedExt(seed_pair_ext, out_file_seed_ext);
  double sxt_dump_time = mytime();
  printElapsed(start_time, sxt_dump_time, "GRASPx::Build info: finished writing seed extensions");
  start_time = mytime();
  // dump read extension
  string out_file_read_ext = workspace_dir + "/" + db_stem + ".rxt";
  db_dump.DumpReadExt(fw_read_ext, re_read_ext, out_file_read_ext);
  double rxt_dump_time = mytime();
  printElapsed(start_time, rxt_dump_time, "GRASPx::Build info: finished writing read extensions");
  start_time = mytime();
  // compute and dump high-scoring k-mers
  unordered_map<string, list<string> > hs_mer;
  db_dump.CreateHighScoreMer(hs_mer);
  string out_file_hs_mer = workspace_dir + "/" + db_stem + ".hsm";
  db_dump.DumpHighScoreMer(hs_mer, out_file_hs_mer);
  double hsmer_time = mytime();
  printElapsed(start_time, hsmer_time, "GRASPx::Build info: finished indexing high-scoring k-mer matches");
  cout << "============================================================" << endl;
  
  
  /*
  // build forward suffix array
  db_seq.BuildSFADefault();
  db_seq.BuildKeyArrayDefault();
  // dump the forward version suffix array
  db_seq.DumpSFA(workspace_dir, db_stem);
  // load in sequence from file for reverse sequences
  
  // create seed-mer map for given database
  DatabaseIndex db_dump(alph_id, seed_len, extd_len);
  
  
  db_dump.BuildSeedmerMap(db_seq);
  
  double sfa_build_time = mytime();
  printElapsed(start_time, sfa_build_time, "build suffix array");
  start_time = mytime();
  
  // create and dump map for reduced sequence to original sequence
  unordered_map<std::string, std::list<std::string> > reduc_alph_map;
  db_dump.CreateReducedMap(reduc_alph_map);
  string out_file = workspace_dir + "/" + db_stem + ".rdm";
  db_dump.DumpReducedMap(reduc_alph_map, out_file);
  
  double rdc_build_time = mytime();
  printElapsed(start_time, rdc_build_time, "build reduced-alphabet mapping");
  start_time = mytime();
  
  // create and dump seed extensions
  unordered_map<string, list<ReadPairType> > seed_pair_ext;
  db_dump.CreateSeedExt(min_seed_coverage, db_seq, db_seq_rev, seed_pair_ext);
  string out_file_seed_ext = workspace_dir + "/" + db_stem + ".sxt";
  db_dump.DumpSeedExt(seed_pair_ext, out_file_seed_ext);
  
  double sde_build_time = mytime();
  printElapsed(start_time, sde_build_time, "build seed extension");
  start_time = mytime();
  
  // create and dump maximal extension reads links
  unordered_map<RIDType, list<OverlapType> > fw_read_ext, re_read_ext;
  db_dump.CreateReadExt(min_ext_coverage, db_seq, db_seq_rev, fw_read_ext, re_read_ext);
  string out_file_read_ext = workspace_dir + "/" + db_stem + ".rxt";
  db_dump.DumpReadExt(fw_read_ext, re_read_ext, out_file_read_ext);
  
  double rde_build_time = mytime();
  printElapsed(start_time, rde_build_time, "build read extension");
  start_time = mytime();
  
  unordered_map<string, list<string> > hs_mer;
  db_dump.CreateHighScoreMer(hs_mer);
  string out_file_hs_mer = workspace_dir + "/" + db_stem + ".hsm";
  db_dump.DumpHighScoreMer(hs_mer, out_file_hs_mer);
  
  double hsm_build_time = mytime();
  printElapsed(start_time, hsm_build_time, "build high-score 3-mer map");
  start_time = mytime();
  */
  return 0;
}
