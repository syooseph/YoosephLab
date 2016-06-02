#include "sequence_build.h"
#include "remap.h"
//#include "database_index.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

using namespace std;

static string workspace_dir;
static string db_file;
static string contig_file;
static string output;
static int seed_size;
static int num_errors;
static double portion_mapped = 0.9;
static bool gap_allowed;

void PrintUsageBuild()  {
  cout << "Usage: graps-map db contig out" << endl;
  cout << "Use \'--help\' for more options" << endl;
  return;
}

bool _cmp_map_reads(MapReadType& m1, MapReadType& m2)  {
  if(m1.rid < m2.rid) {
    return true;
  }
  return false;
}

int main(int argc, char** argv)  {
  // reading options
  boost::program_options::options_description desc("List of options");
  desc.add_options()
      ("help", "print the help message")
      ("db", boost::program_options::value<string>(&db_file), "short-peptide reads (in FASTA format)")
      ("contig", boost::program_options::value<string>(&contig_file), "contigs output by grasp-assemble")
      ("out", boost::program_options::value<string>(&output), "read headers mapped to the contigs")
      ("work_space", boost::program_options::value<string>(&workspace_dir)->default_value("WorkSpace"), "working directory for indexing file dump")
      ("seed_size", boost::program_options::value<int>(&seed_size)->default_value(10), "seed size for initiating local alignment")
      ("num_errors", boost::program_options::value<int>(&num_errors)->default_value(2), "num errors allowed for mapping")
      //("gap_allowed", boost::program_options::value<bool>(&gap_allowed)->default_value(false), "whether allow gap")
      ("portion_mapped", boost::program_options::value<double>(&portion_mapped), "minimum portion of read needs to be aligned (default 0.9)")
  ;
  boost::program_options::positional_options_description pos_opt;
  pos_opt.add("db", 1);
  pos_opt.add("contig", 1);
  pos_opt.add("out", 1);
  boost::program_options::variables_map vm;
  boost::program_options::store(
      boost::program_options::command_line_parser(argc, argv).
      options(desc).positional(pos_opt).run(), vm
  );
  boost::program_options::notify(vm);
  if(vm.count("help"))  {
    cout << "Usage: graps-map db contig out" << endl;
    cout << desc << endl; 
    return 0;
  }
  // check options validity
  if(!boost::filesystem::exists(db_file))  {
    cout << db_file << endl;
    cout << "Error: grasp-map: db_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::exists(contig_file))  {
    cout << contig_file << endl;
    cout << "Error: grasp-map: contig_file does not exist." << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(!boost::filesystem::exists(workspace_dir))  {
    cout << "Error: grasp-map: working space does not exist. " << workspace_dir << endl;
    cout << "Please use \'--help\' for more details." << endl;
    exit(0);
  }
  if(num_errors < 0 || num_errors > 3)  {
    cout << "Error: grasp-map: num_errors out of range (allowed range: 0-3)." << endl;
    exit(0);
  }
  if(portion_mapped < 0.6 || portion_mapped > 1.0)  {
    cout << "Error: grasp-map: portion_mapped out of range (allowed range: 0.6-1.0)." << endl;
    exit(0);
  }
  if(seed_size < 5 || portion_mapped > 30)  {
    cout << "Error: grasp-map: seed_size out of range (allowed range: 5-30)." << endl;
    exit(0);
  }
  // load the sequence and suffix array
  SequenceBuild db_seq(db_file);
  string db_stem = db_seq.GetFileStem(db_file);
  db_seq.LoadSFA(workspace_dir, db_stem);
  
  // DEBUG BEGIN FOR NEW SUFFIX ARRAY
  /*
  DatabaseIndex kmer_loader;
  unordered_map<string, list<ReadPairType> > ext_seed;
  string idx_file = "/usr/local/depot/projects/SPA/czhong/GRASPx_current/WorkSpace/stool.faa.sxt";
  kmer_loader.LoadSeedExt(idx_file, ext_seed);
  for(auto it = ext_seed.begin(); it != ext_seed.end(); ++ it) {
    string mer = it->first;
    pair<IdxType, IdxType> bound = db_seq.SearchSFA(mer);
    IdxType coverage = bound.second - bound.first + 1;
    cout << (unsigned long int) coverage << endl;
  }
  return 0;
  */
  // DEBUG END  
  //cout << "Finish loading reads" << endl;
  // reads in the contig sequences
  ReMap re_aln;
  list<ContigType> loaded_contigs;
  re_aln.LoadContigs(contig_file, loaded_contigs);
  //cout << "Finish loading contigs" << endl;
  // remap the reads back to the contigs
  list<MapReadType> identified_reads;
  re_aln.MapToContig(db_seq, loaded_contigs, seed_size, portion_mapped, num_errors, gap_allowed, identified_reads);
  identified_reads.sort(_cmp_map_reads);
  // output the headers
  bool write_to_file = false;
  ofstream out_fh;
  if(vm.count("out"))  {
    out_fh.open(output.c_str(), ios_base::out);
    if(out_fh.good())  {
      write_to_file = true;
    }
  } 
  if(write_to_file)  {
    for(auto it = identified_reads.begin(); it != identified_reads.end(); ++ it) {
      out_fh << (int) (it->rid) << " " << db_seq.GetHeader(it->rid) << "  " << re_aln.GetContigInfo(it->contig_id) << endl;
    }
    out_fh.close();
  } else  {
    for(auto it = identified_reads.begin(); it != identified_reads.end(); ++ it) {
      cout << (int) (it->rid) << " " << db_seq.GetHeader(it->rid) << "  " << re_aln.GetContigInfo(it->contig_id) << endl;
    } 
  }
  return 0;
}
